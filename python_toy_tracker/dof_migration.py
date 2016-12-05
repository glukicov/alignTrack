import toytraceback as ttb
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys
import math
import scipy.stats as stats
import os.path
import getopt
import ROOT

track_count = 0
fit_count = 0
fits_attempted = 10000
chi_sq_cut = 2
missing_straws = False
input_dir = ""
output_dir = ""


argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "h:i:t:c:o:m:f:", ["help", "input_dir=", "track_count=", "chi_sq_cut=", "output_dir=", "missing_straws=", "fits_attempted="])
except getopt.GetoptError:
    print "tree_reader.py -i <input_dir> -t <track_count> -f <fit_count> -c <chi_sq_cut> -o <output_dir> -m <missing_straws>"
    sys.exit(2)

# Set parameters according to arguments
for opt, arg in opts:

    if opt in ("-h", "--help"):
        print "tree_reader.py -i <input_dir> -t <track_count> -c <chi_sq_cut> -o <output_dir>"
        sys.exit()
    elif opt in ("-i", "--input_dir"):
        input_dir = arg
    elif opt in ("-t", "--track_count"):
        track_count = int(arg)
    elif opt in ("-c", "--chi_sq_cut"):
        chi_sq_cut = float(arg)
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    elif opt in ("-m", "--missing_straws"):
        if arg == "True":
            missing_straws = True
        elif arg == "False":
            missing_straws = False
    elif opt in ("-f", "--fits_attempted"):
        fits_attempted = int(arg)

# String for input file name
straw_status = "missed_straws" if missing_straws else "all_straws"
input_file = input_dir + "smear_" + straw_status + "_" + str(track_count) + "_track_" + str(fits_attempted) + "_fit.root"

print input_file

# Get input file
f = ROOT.TFile.Open(input_file, "read")

# Get tree of straw hit events
tree = f.Get("event_tree")

# Get number of fits in tree, by examining last entry in tree 
entry_count = tree.GetEntries()
tree.GetEntry(entry_count - 1)
fit_count = tree.fitNum + 1

# List for number of degrees of freedom in fits
fit_dof = []

print "Total Fits:", str(fit_count)
        
for i in xrange(fit_count):
    
    # Lists of hit distances for fit
    true_hit_distances = []
    fitted_hit_distances = []

    if (i % 100) == 0:
        print i

    # Get entries with correct fit number in tree
    fit_cut_tree = tree.CopyTree('fitNum==' + str(i))

    # Get number of entries for this fit, and loop over them
    fit_cut_tree_entry_count = fit_cut_tree.GetEntries() 
            
    for j in xrange(fit_cut_tree_entry_count):
        
        fit_cut_tree.GetEntry(j)

        # Append hit distances for this entry
        true_hit_distances.append(fit_cut_tree.trueHitDistance)
        fitted_hit_distances.append(fit_cut_tree.fittedHitDistance)


    # Create array of hit distance uncertainties, taken to be equal to hit resolution of detector    
    sigmas = np.array([ttb.hit_resolution for n in xrange(fit_cut_tree_entry_count)])

    # Calculate reduced value of chi-squared from hit distances
    chi_sq_red = np.sum(((np.array(true_hit_distances) - np.array(fitted_hit_distances)) / sigmas)**2) / (fit_cut_tree_entry_count - ((2 * track_count) + 1)) 

    # Check if chi-squared below cut, then record number of degrees of freedom for fit
    if (chi_sq_red < chi_sq_cut):
        fit_dof.append(fit_cut_tree_entry_count - ((2 * track_count) + 1))

    fit_cut_tree.IsA().Destructor(fit_cut_tree) # To avoid stupid bloody memory leak



# Calculate stats for list of degrees of freedoms
dof_mean = np.mean(fit_dof)
dof_std_dev = np.std(fit_dof)
dof_data_count = len(fit_dof)

# Further stats for histogramming
dof_max = max(fit_dof)
dof_min = min(fit_dof)
dof_range = dof_max - dof_min

# Plot histogram of alignment errors
plt.hist(fit_dof, bins=dof_range + 3, range=[dof_min - 1, dof_max + 2], color='blue')
plt.title("Plot of Module Misalignment Distances for Fits With $\chi^2_{red} < " + str(chi_sq_cut) + "$, With" + str(track_count) + " Tracks")
plt.xlabel("Number of Degrees of Freedom in Fit")
plt.ylabel("Frequency")

# Entries to put in legend, with relevent statistical quantities
number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(dof_data_count)))
mean_patch = mpatches.Patch(color='none', label=("$\mu$: " + str(dof_mean)))
err_mean_patch = mpatches.Patch(color='none', label=("$SE(\mu)$: " + str(dof_std_dev / math.sqrt(dof_data_count))))
std_patch = mpatches.Patch(color='none', label=("$\sigma$: " + str(dof_std_dev)))
err_std_patch = mpatches.Patch(color='none', label=("$SE(\sigma)$: " + str(dof_std_dev / math.sqrt((2 * dof_data_count) - 2))))

# Get legend handles and labels, and plot legend
dof_handles = [number_patch, mean_patch, err_mean_patch, std_patch, err_std_patch]
dof_labels = [number_patch.get_label(), mean_patch.get_label(), err_mean_patch.get_label(), std_patch.get_label(), err_std_patch.get_label()]
plt.legend(dof_handles, dof_labels, loc='best', frameon=False)

# Save plot
plt.savefig(output_dir + straw_status + "_" + str(track_count) + "t_dof.pdf")
plt.clf()
