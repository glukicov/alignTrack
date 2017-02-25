# 
# Script to load fitted, true parameters from Root file, then plot the difference
# between these values, normalised by the fitted parameter uncertainty, against
# the label number of these parameters.
#
# John Smeaton 07/02/2017
#  

import ROOT
import sys
import getopt
import os
import numpy as np
import matplotlib.pyplot as plt
import millepede_utils


# Get system arguments, define string showing help message
argv = sys.argv[1:]
helpstring = "compareLabelsRoot.py -i <input_file> -o <output_dir>"

input_filename = "" # Define filename string

# Get options, arguments
try:
    opts, args = getopt.getopt(argv, "h:i:", ["help", "input_file"])
except getopt.GetoptError:
    print helpstring
    sys.exit(2)

# Set input filename, if given in console argument
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print helpstring
        sys.exit(2)
    elif opt in ("-i", "--input_file"):
        input_filename = os.path.abspath(arg)


# Open file, and get TTree
f = ROOT.TFile.Open(input_filename, "read")
tree = f.Get("paramTree")

# Get subsets of this tree, with true values (fitType==0), and fitted values 
# using the inversion method for fitting (fitType==1).
fit_tree = tree.CopyTree('fitType==' + str(1.0))
true_tree = tree.CopyTree('fitType==' + str(0.0))

# Dictionaries for true and fitted parameter values, and fitted parameter errors,
# with parameter labels acting as keys.
param_true_vals = dict()
param_fit_vals = dict()
param_errors = dict()

# Loop across tree entries
for i in xrange(fit_tree.GetEntries()):

    # Get true parameter value, and add to dictionary
    true_tree.GetEntry(i)
    param_true_vals[true_tree.label] = true_tree.paramValue

    # Get fitted parameter value, error, and add to dictionaries
    fit_tree.GetEntry(i)
    param_fit_vals[fit_tree.label] = fit_tree.paramValue
    param_errors[fit_tree.label] = fit_tree.paramError


# Lists of normalised parameter differences for plane displacements, and for 
# velocity deviations
fit_error_plane_y_disp = []
fit_error_vel_dev = []

# Lists for parameter labels, matched to differences stored above
label_plane_y_disp = []
label_vel_dev = []

# Lists for true parameter values
true_plane_y_disp = []
true_vel_dev = []

# Lists for fitted parameter values
fit_plane_y_disp = []
fit_vel_dev = []



# Loop across all recorded parameter labels
for label in param_true_vals.iterkeys():
    
    # Check label value, then calculate normalised parameter difference, adding 
    # values to appropriate lists
    if label < 500:
        label_plane_y_disp.append(label)
        true_plane_y_disp.append(param_true_vals[label])
        fit_plane_y_disp.append(param_fit_vals[label])
        fit_error_plane_y_disp.append((param_fit_vals[label] - param_true_vals[label]) / param_errors[label])
    else:
        label_vel_dev.append(label)
        true_vel_dev.append(param_true_vals[label])
        fit_vel_dev.append(param_fit_vals[label])
        fit_error_vel_dev.append((param_fit_vals[label] - param_true_vals[label]) / param_errors[label])



# Plot scatter graph of plane displacements against parameter label
plt.plot(label_plane_y_disp, true_plane_y_disp, 'b.')
plt.title("True Plane Y-Displacement Parameter Values")
plt.ylabel("Parameter Value")
plt.xlabel("Parameter Label")
plt.show()

# Plot scatter graph of plane displacements against parameter label
plt.plot(label_vel_dev, true_vel_dev, 'b.')
plt.title("True Velocity Deviation Parameter Values")
plt.ylabel("Parameter Value")
plt.xlabel("Parameter Label")
plt.show()


# Plot scatter graph of plane displacements against parameter label
plt.plot(label_plane_y_disp, fit_plane_y_disp, 'b.')
plt.title("Fitted Plane Y-Displacement Parameter Values")
plt.ylabel("Parameter Value")
plt.xlabel("Parameter Label")
plt.show()

# Plot scatter graph of plane displacements against parameter label
plt.plot(label_vel_dev, fit_vel_dev, 'b.')
plt.title("Fitted Velocity Deviation Parameter Values")
plt.ylabel("Parameter Value")
plt.xlabel("Parameter Label")
plt.show()


# Plot scatter graph of plane displacements against parameter label
plt.plot(label_plane_y_disp, np.array(fit_plane_y_disp) - np.array(true_plane_y_disp), 'b.')
plt.title("Unnormalised Differences Between Fitted, True \n Plane Y-Displacement Parameter Values")
plt.ylabel("Parameter Value")
plt.xlabel("Parameter Label")
plt.show()

# Plot scatter graph of plane displacements against parameter label
plt.plot(label_vel_dev, np.array(fit_vel_dev) - np.array(true_vel_dev), 'b.')
plt.title("Unnormalised Differences Between Fitted, True \n Velocity Deviation Parameter Values")
plt.ylabel("Parameter Value")
plt.xlabel("Parameter Label")
plt.show()


# Plot scatter graph of plane displacement differences against parameter label
plt.plot(label_plane_y_disp, fit_error_plane_y_disp, 'b.')
plt.title("Differences Between Fitted, True Plane Y-Displacement Parameters, \n Divided by Fitted Parameter Uncertainty")
plt.ylabel("Normalised Difference")
plt.xlabel("Parameter Label")
plt.show()

# Plot scatter graph of velocity deviation differences against parameter label
plt.plot(label_vel_dev, fit_error_vel_dev, 'b.')
plt.title("Differences Between Fitted, True Velocity Deviation Parameters, \n Divided by Fitted Parameter Uncertainty")
plt.ylabel("Normalised Difference")
plt.xlabel("Parameter Label")
plt.show()
