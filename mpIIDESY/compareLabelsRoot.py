# 
# Script to load fitted, true parameters from Root file, then plot the difference
# between these values, normalised by the fitted parameter uncertainty, against
# the label number of these parameters.
#
# John Smeaton 20/03/2017
#  

import ROOT
import sys
import getopt
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib.pylab as pylab
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
fit_error_plane_rot = []

# Lists for parameter labels, matched to differences stored above
label_plane_y_disp = []
label_vel_dev = []
label_plane_rot = []

# Lists for true parameter values
true_plane_y_disp = []
true_vel_dev = []
true_plane_rot = []

# Lists for fitted parameter values
fit_plane_y_disp = []
fit_vel_dev = []
fit_plane_rot = []


# Loop across all recorded parameter labels
for label in param_true_vals.iterkeys():
    
    # Check label value, then calculate normalised parameter difference, adding 
    # values to appropriate lists
    if label < 500:
        label_plane_y_disp.append(label)
        true_plane_y_disp.append(param_true_vals[label])
        fit_plane_y_disp.append(param_fit_vals[label])
        fit_error_plane_y_disp.append(param_errors[label])
    elif label < 1000:
        label_vel_dev.append(label)
        true_vel_dev.append(param_true_vals[label])
        fit_vel_dev.append(param_fit_vals[label])
        fit_error_vel_dev.append(param_errors[label])
    else:
        label_plane_rot.append(label)
        true_plane_rot.append(param_true_vals[label])
        fit_plane_rot.append(param_fit_vals[label])
        fit_error_plane_rot.append(param_errors[label])


plane_y_disp_chisq = np.sum(((np.array(fit_plane_y_disp) - np.array(true_plane_y_disp)) / np.array(fit_error_plane_y_disp))**2) / len(true_plane_y_disp)
plane_rot_chisq = np.sum(((np.array(fit_plane_rot) - np.array(true_plane_rot)) / np.array(fit_error_plane_rot))**2) / len(true_plane_rot)
vel_dev_chisq = np.sum(((np.array(fit_vel_dev) - np.array(true_vel_dev)) / np.array(fit_error_vel_dev))**2) / len(true_vel_dev)

plane_y_disp_chisq_patch = mpatch.Patch(color='none', label=("Param $\\chi^2_{red}$: " + "{0:.2f}".format(plane_y_disp_chisq)))
plane_rot_chisq_patch = mpatch.Patch(color='none', label=("Param $\\chi^2_{red}$: " + "{0:.2f}".format(plane_rot_chisq)))
vel_dev_chisq_patch = mpatch.Patch(color='none', label=("Param $\\chi^2_{red}$: " + "{0:.2f}".format(vel_dev_chisq)))


textparams = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(textparams)

# # Plot scatter graph of true plane displacements against parameter label
# plt.plot(label_plane_y_disp, true_plane_y_disp, 'b.')
# plt.title("True Plane Y-Displacement Parameter Values")
# plt.ylabel("Parameter Value")
# plt.xlabel("Parameter Label")
# plt.show()

# # Plot scatter graph of true rotations against parameter label
# plt.plot(label_plane_rot, true_plane_rot, 'b.')
# plt.title("True Plane Rotation Parameter Values")
# plt.ylabel("Parameter Value")
# plt.xlabel("Parameter Label")
# plt.show()

# # Plot scatter graph of true drift velocities against parameter label
# plt.plot(label_vel_dev, true_vel_dev, 'b.')
# plt.title("True Velocity Deviation Parameter Values")
# plt.ylabel("Parameter Value")
# plt.xlabel("Parameter Label")
# plt.show()


# # Plot scatter graph of fitted plane displacements against parameter label
# plt.errorbar(label_plane_y_disp, fit_plane_y_disp, xerr=0, yerr=fit_error_plane_y_disp, color='b', fmt='.')
# plt.title("Fitted Plane Y-Displacement Parameter Values")
# plt.ylabel("Parameter Value")
# plt.xlabel("Parameter Label")
# plt.show()

# # Plot scatter graph of fitted plane displacements against parameter label
# plt.errorbar(label_plane_rot, fit_plane_rot, xerr=0, yerr=fit_error_plane_rot, color='b', fmt='.')
# plt.title("Fitted Plane Rotation Parameter Values")
# plt.ylabel("Parameter Value")
# plt.xlabel("Parameter Label")
# plt.show()

# # Plot scatter graph of fitted velocity deviations against parameter label
# plt.errorbar(label_vel_dev, fit_vel_dev, xerr=0, yerr=fit_error_vel_dev, color='b', fmt='.')
# plt.title("Fitted Velocity Deviation Parameter Values")
# plt.ylabel("Parameter Value")
# plt.xlabel("Parameter Label")
# plt.show()


# # Plot scatter graph of fitted plane displacements against parameter label
# plt.errorbar(true_plane_y_disp, fit_plane_y_disp, xerr=0, yerr=fit_error_plane_y_disp, color='b', fmt='.')
# plt.title("Fitted vs. True Plane Y-Displacement Parameter Values")
# plt.ylabel("Fitted Y-Displacement")
# plt.xlabel("True Y-Displacement")
# plt.show()

# # Plot scatter graph of fitted plane displacements against parameter label
# plt.errorbar(true_plane_rot, fit_plane_rot, xerr=0, yerr=fit_error_plane_rot, color='b', fmt='.')
# plt.title("Fitted vs. True Plane Rotation Parameter Values")
# plt.ylabel("Fitted Plane Rotation")
# plt.xlabel("True Plane Rotation")
# plt.show()

# # Plot scatter graph of fitted velocity deviations against parameter label
# plt.errorbar(true_vel_dev, fit_vel_dev, xerr=0, yerr=fit_error_vel_dev, color='b', fmt='.')
# plt.title("Fitted vs. True Velocity Deviation Parameter Values")
# plt.ylabel("Fitted Velocity Deviation")
# plt.xlabel("True Velocity Deviation")
# plt.show()


# Plot scatter graph of unnormalised plane displacement fitted - true values 
plt.errorbar(label_plane_y_disp, np.array(fit_plane_y_disp) - np.array(true_plane_y_disp), xerr=0, yerr=fit_error_plane_y_disp, color='b', fmt='.')
plt.title("Differences Between Fitted, True \n Plane Y-Displacements ")
plt.ylabel("(Fitted - True) Plane Displacement")
plt.xlabel("Parameter Label")
plt.legend(handles=[plane_y_disp_chisq_patch], loc='best')
plt.grid()
plt.show()

# Plot scatter graph of unnormalised plane displacement fitted - true values 
plt.errorbar(label_plane_rot, np.array(fit_plane_rot) - np.array(true_plane_rot), xerr=0, yerr=fit_error_plane_rot, color='b', fmt='.')
plt.title("Differences Between Fitted, True \n Plane Rotations")
plt.ylabel("(Fitted - True) Plane Rotation")
plt.xlabel("Parameter Label")
plt.legend(handles=[plane_rot_chisq_patch], loc='best')
plt.grid()
plt.show()

# Plot scatter graph of unnormalised velocity deviations fitted - true values 
plt.errorbar(label_vel_dev, np.array(fit_vel_dev) - np.array(true_vel_dev), xerr=0, yerr=fit_error_vel_dev, color='b', fmt='.')
plt.title("Differences Between Fitted, True \n Velocity Deviations")
plt.ylabel("(Fitted - True) Velocity Deviation")
plt.xlabel("Parameter Label")
plt.legend(handles=[vel_dev_chisq_patch], loc='best')
plt.grid()
plt.show()


# # Plot scatter graph of normalised plane displacement fitted - true values 
# plt.plot(label_plane_y_disp, (np.array(fit_plane_y_disp) - np.array(true_plane_y_disp)) / np.array(fit_error_plane_y_disp), 'b.')
# plt.title("Differences Between Fitted, True Plane Y-Displacement Parameters, \n Divided by Fitted Parameter Uncertainty")
# plt.ylabel("(Fitted - True) Plane Displacement - Unc. Normalised")
# plt.xlabel("Parameter Label")
# plt.show()

# # Plot scatter graph of normalised plane displacement fitted - true values 
# plt.plot(label_plane_rot, (np.array(fit_plane_rot) - np.array(true_plane_rot)) / np.array(fit_error_plane_rot), 'b.')
# plt.title("Differences Between Fitted, True Plane Rotation Parameters, \n Divided by Fitted Parameter Uncertainty")
# plt.ylabel("(Fitted - True) Plane Rotation - Unc. Normalised")
# plt.xlabel("Parameter Label")
# plt.show()


# # Plot scatter graph of normalised velocity deviations fitted - true values 
# plt.plot(label_vel_dev, (np.array(fit_vel_dev) - np.array(true_vel_dev)) / np.array(fit_error_vel_dev), 'b.')
# plt.title("Differences Between Fitted, True Velocity Deviation Parameters, \n Divided by Fitted Parameter Uncertainty")
# plt.ylabel("(Fitted - True) Velocity Deviation - Unc. Normalised")
# plt.xlabel("Parameter Label")
# plt.show()


# # Plot scatter graph of fitted plane displacements against parameter label
# plt.errorbar(true_plane_rot, np.array(fit_plane_rot) - np.array(true_plane_rot), xerr=0, yerr=fit_error_plane_rot, color='b', fmt='.')
# plt.grid()
# plt.title("Difference Between Fitted, True Plane Rotations, \n As a Function of True Plane Rotations")
# plt.ylabel("(Fitted - True) Plane Rotation / rad")
# plt.xlabel("True Plane Rotation / rad")
# plt.legend(handles=[plane_rot_chisq_patch], loc='best')
# plt.show()
