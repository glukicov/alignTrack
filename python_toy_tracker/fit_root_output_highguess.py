import toytraceback as ttb
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import random
import sys
import math
import scipy.stats as stats
import time
import os.path
import getopt
from scipy.optimize import curve_fit 

from rootpy.tree import Tree, TreeModel
from rootpy.tree import IntCol, DoubleCol, DoubleArrayCol, ObjectCol
from rootpy.io import root_open
from rootpy import stl
import rootpy.compiled

# Script to carry out multiple fits with a specified number of tracks, recording results in a ROOT TTree.

# Model for each entry in tree
class HitEvent(TreeModel):
    
    # Integers for number of fit, and number of event
    fitNum = IntCol()
    eventNum = IntCol()
    
    # Doubles for fitted and true alignment for this fit
    fittedAlignment = DoubleCol()
    trueAlignment = DoubleCol()
    
    # Doubles for track parameters
    fittedTrackGrad = DoubleCol()
    fittedTrackInt = DoubleCol()
    trueTrackGrad = DoubleCol()
    trueTrackInt = DoubleCol()
    
    # Doubles for true and fitted hit distances
    trueHitDistance = DoubleCol()
    fittedHitDistance = DoubleCol()

    # Indexes for position of wire hit
    moduleNum = IntCol()
    planeNum = IntCol()
    layerNum = IntCol()
    wireNum = IntCol()

# Get start CPU time, and wall time (in seconds), to ensure script doesn't overrun batch system limits
start_wall_time = time.time()
start_cpu_time = time.clock()

track_count = 5 # Number of tracks to fit
module_alignment = 1.5 # Alignment of first module (in mm)
fit_count = 10 # Number of times to run fit

# Whether to include smearing of residuals, and finite straw size
detector_smear = True
detector_missing_hits = False

# Time limits for script on batch system (in s), with time buffer to ensure partial results are processed. 
wall_time_limit = 345600
cpu_time_limit = 259200
time_buffer = 100

# Directory where plots are to be saved
directory = os.path.expanduser("~/")

# Dictionary of CPU and wall times for different UCL Batch queues
queue_time_dict = {"short": [3600, 3600], "medium": [28800, 86400], "long": [259200, 345600]}

# Get arguments and options from console

print str(sys.argv)

argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "h:f:t:s:m:q:a:o:", ["help", "fit_count=", "track_count=", "smearing=", "missing_hits=", "batch_queue=", "module_alignment=", "output_dir="])
except  getopt.GetoptError:
    print "multi_fit.py -f <fit_count> -t <track_count> -s <smearing> -m <missing_hits> -q <batch_queue> -a <module_alignment> -o <output_dir>"
    sys.exit(2)

# Set parameters according to arguments
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print "multi_fit.py -f <fit_count> -t <track_count> -s <smearing> -m <missing_hits> -q <batch_queue> -a <module_alignment> -o <output_dir>"
        sys.exit()
    elif opt in ("-f", "--fit_count"):
        fit_count = int(arg)
    elif opt in ("-t", "--track_count"):
        track_count = int(arg)
    elif opt in ("-s", "--smearing"):
        if (arg == "True"):
            detector_smear = True
        elif (arg == "False"):
            detector_smear = False
        else:
            print "-s, --smearing: Must use 'True' or 'False'"
            sys.exit()
    elif opt in ("-m", "--missing_hits"):
        if (arg == "True"):
            detector_missing_hits = True
        elif (arg == "False"):
            detector_missing_hits = False
        else:
            print "-m, --missing_hits: Must use 'True' or 'False'"
            sys.exit()
    elif opt in ("-q", "--batch_queue"):
        if (arg in ("short", "medium", "long")):
            cpu_time_limit = queue_time_dict[arg][0]
            wall_time_limit = queue_time_dict[arg][1]
        else:
            print "-q, --batch_queue: Must use 'short', 'medium', or 'long'"
            sys.exit()
    elif opt in ("-a", "-module_alignment"):
        module_alignment = float(arg)
    elif opt in ("-o", "-output_dir"):
        directory = arg
    
# Set file name for output file, then open root file and create tree 
smear_string = "smear" if detector_smear else "no_smear"
straw_string = "missed_straws" if detector_missing_hits else "all_straws"
file_name = directory + smear_string + "_" + straw_string + "_" + str(track_count) + "_track_" + str(fit_count) + "_fit.root"

print file_name
f = root_open(file_name, "recreate")
t = Tree("event_tree", model=HitEvent)

            
print "Carrying out", fit_count, "fits, with", track_count, "tracks:"  
        
# Carry out requested number of fits
for k in xrange(fit_count):
    
    print k
    
    # FITTING

    # Create instance of detector, and set alignment to desired value
    true_detector = ttb.Detector(smearing=detector_smear, finite_straws=detector_missing_hits)
    true_detector.set_module_x_align(1, module_alignment)

    # Coordinates at beginning, end of track
    x_bottom = [random.uniform(ttb.track_boundary_x_lower, ttb.track_boundary_x_upper) for i in xrange(track_count)]
    x_top = [random.uniform(ttb.track_boundary_x_lower, ttb.track_boundary_x_upper) for i in xrange(track_count)]
    y_bottom = [ttb.track_boundary_y_lower for i in xrange(track_count)]
    y_top = [ttb.track_boundary_y_upper for i in xrange(track_count)]
    
        
    # Calculate track gradient and x-intercept from track coordinates
    gradient = [((x_top[i] - x_bottom[i]) / (y_top[i] - y_bottom[i])) for i in xrange(track_count)]
    intercept = [x_bottom[i] - (gradient[i] * y_bottom[i]) for i in xrange(track_count)]
        
    # Set up track, and get pos of beginning, end points
    true_track = [ttb.Track(gradient[i], intercept[i]) for i in xrange(track_count)]
        
    # Get wire hits in each layer, then assign these to detector
    events = []
    for j in xrange(track_count): 
        wire_hits = ttb.closest_hit_wires(true_detector, true_track[j])
        event = ttb.DetectorHitEvent(true_track[j], wire_hits, j)
        events.append(event)
            

    # New detector object for fitting
    fitting_detector = ttb.Detector(smearing=detector_smear, finite_straws=detector_missing_hits)
    fitting_detector.set_events(events)

    # Get x-displacement of all closest approached wires
    hit_rads = []
    for event in events:
        for wire_hit in event.wire_hits:
            hit_rads.append(wire_hit.get_hit_dist())
                
    # Get list of keys for all struck wires, by all tracks
    wire_keys = fitting_detector.get_hit_keys()

    guess = [0.0 for i in xrange((2 * track_count) + 1)] # Initial guess for fitted track gradient and intercept, and module alignment (spurious convergence without this)
    guess[-1] = 3.0
    

    fit_sigmas = [ttb.hit_resolution]*len(hit_rads) # Uncertainties on hit radii, from detector hit resolution

    # Finds values for track gradient and intercept, by fitting to wire hit x-displacements
    popt, pcov = curve_fit(fitting_detector.get_hit_radius, wire_keys, hit_rads, p0=guess, sigma=fit_sigmas)


    # Loop across events in this fit, creating one entry in tree for each event
    for event_num in xrange(len(events)):


        # Get wire keys for this event only
        filtered_wire_keys = [wire_key for wire_key in wire_keys if int(wire_key.split('-')[0]) == event_num]
        
        # Get arrays of true and fitted hit radii for this event
        true_hit_rads = [wire_hit.get_hit_dist() for wire_hit in events[event_num].wire_hits]
        fitted_hit_rads = fitting_detector.get_hit_radius(filtered_wire_keys, popt)

        # Records hit distances to arrays, recording zero if straw was not hit.
        for j in xrange(len(true_hit_rads)):
            # Record event number and fit number for entry in tree
            t.eventNum = event_num
            t.fitNum = k

            # Record true fitted alignment
            t.trueAlignment = module_alignment
            t.fittedAlignment = popt[-1]

            # Record true and fitted track parameters
            t.trueTrackGrad = events[event_num].track.get_gradient()
            t.trueTrackInt = events[event_num].track.get_intercept()
            t.fittedTrackGrad = popt[2 * event_num]
            t.fittedTrackInt = popt[(2 * event_num) + 1] 

            # Record true and fitted hit distances
            t.trueHitDistance = true_hit_rads[j]
            t.fittedHitDistance = fitted_hit_rads[j]
            
            # Record numbers indexing position of struck wire
            t.moduleNum = events[event_num].wire_hits[j].module_num
            t.planeNum = events[event_num].wire_hits[j].plane_num
            t.layerNum = events[event_num].wire_hits[j].layer_num
            t.wireNum =  events[event_num].wire_hits[j].wire_num

            # Fill tree
            t.Fill()
 
    # Stop carrying out fits, if close to time limit
    if (time.time() - start_wall_time) > (wall_time_limit - time_buffer):
        fit_count = k + 1
        break
    if (time.clock() - start_cpu_time) > (cpu_time_limit - time_buffer):
        fit_count = k + 1
        break


# Write tree to file, then close
f.Write()
f.Close()

print ""
print("Closed TFile")
print ""

# Get end time of script, then output    
end_wall_time = time.time()
end_cpu_time = time.clock()
print "Completed in " + str(end_wall_time - start_wall_time) + "s."
print "CPU Time Elapsed: " + str(end_cpu_time - start_cpu_time) + "s."
print ""
