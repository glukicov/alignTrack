import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab
import numpy as np
import sys
import getopt
import os


# Filenames for txt outputs of true, fitted parameters
true_params_filename = ""
fitted_params_filename = ""

# Get console arguments
argv = sys.argv[1:]

# Set up options for arguments
helpstring = "compareParams.py -f <fitted_params> -t <true_params>"
try:
    opts, args = getopt.getopt(argv, "h:t:f:", ["help", "fitted_params", "true_params"])
except getopt.GetoptError:
    print helpstring
    sys.exit(2)

# Get fitted, true parameter filenames from arguments
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print helpstring
        sys.exit()
    elif opt in ("-f", "-fitted_params"):
        fitted_params_filename = arg
    elif opt in ("-t", "-true_params"):
        true_params_filename = arg

# Dictionaries to contain values of parameters and errors. Uses parameter labels as keys
true_params = dict()
fitted_params = dict()
fitted_param_errors = dict()


# Open text file of true parameters, iterating across all lines
with open(true_params_filename, 'r') as true_f:
    for line in true_f.readlines():

        # Get items in line, then get parameter value, if first item in line is a label
        items = line.split()
        try: 
            int(items[0])
            true_params[int(items[0])] = float(items[1])
        except Exception:
            pass

# Open text file of true parameters, iterating across all lines
with open(fitted_params_filename, 'r') as fitted_f:
    for line in fitted_f.readlines():
       
        # Get items in line, then get parameter and error values, if first item in line is a label
        items = line.split()        
        try: 
            int(items[0])
            fitted_params[int(items[0])] = float(items[1])
            fitted_param_errors[int(items[0])] = float(items[4])
        except Exception:
            pass
    
# Lists for parameter values of plane displacement and drift velocity deviation
true_drift_vel = []
true_vert_dev = []

# Lists for difference between fitted and true parameter values, divided by uncertainties in fitted parameter uncertainty
unc_norm_errors_drift_vel = []
unc_norm_errors_vert_dev = []

# Lists for difference between fitted and true parameter values, divided by uncertainties in fitted parameter uncertainty
val_norm_errors_drift_vel = []
val_norm_errors_vert_dev = []

# Absolute differences between fitted, true parameter values
abs_errors_vel = []
abs_errors_dev = []


# Print column headers
if (len(fitted_param_errors) > 0):
    print ""
    print "{:<10} {:<15} {:<15} {:<15} {:<22} {:<22}".format("Label", "True Param.", "Fitted Param.", "Fit Error", "Unc. Norm. Difference", "Val. Norm. Difference")
    print ""
else:
    print ""
    print "{:<10} {:<15} {:<15} {:<22}".format("Label", "True Param.", "Fitted Param.", "Val. Norm. Difference")
    print ""
    


# Iterate across all keys in dictionary
for label in sorted(true_params.iterkeys()):

    # Print blank line between PDD, DVD outputs 
    if label == 501: print ""

    if (len(fitted_param_errors) > 0):

        # Print parameter label, true and fitted values, fitted uncertainty, and normalised difference between fitted, true values, padding appropriately
        print "{:<10} {:< 15} {:< 15} {:< 15} {:< 22} {:< 22}".format(label, true_params[label], fitted_params[label], fitted_param_errors[label], (fitted_params[label] - true_params[label]) / fitted_param_errors[label], (fitted_params[label] - true_params[label]) / fitted_params[label])

        # Add normalised difference, true parameter values to appropriate lists, checking label value to see if entry corrsponds to displacement of velocity deviation
        if label < 500:
            unc_norm_errors_vert_dev.append((fitted_params[label] - true_params[label]) / fitted_param_errors[label])
            val_norm_errors_vert_dev.append((fitted_params[label] - true_params[label]) / fitted_params[label])
            abs_errors_dev.append(fitted_params[label] - true_params[label])
            true_vert_dev.append(true_params[label])  
        else:
            unc_norm_errors_drift_vel.append((fitted_params[label] - true_params[label]) / fitted_param_errors[label])
            val_norm_errors_drift_vel.append((fitted_params[label] - true_params[label]) / fitted_params[label])
            abs_errors_vel.append(fitted_params[label] - true_params[label])
            true_drift_vel.append(true_params[label])


    else:

        # Print parameter label, true and fitted values, fitted uncertainty, and normalised difference between fitted, true values, padding appropriately
        print "{:<10} {:< 15} {:< 15} {:< 22}".format(label, true_params[label], fitted_params[label], (fitted_params[label] - true_params[label]) / fitted_params[label])

        # Add normalised difference, true parameter values to appropriate lists, checking label value to see if entry corrsponds to displacement of velocity deviation
        if label < 500:
            val_norm_errors_vert_dev.append((fitted_params[label] - true_params[label]) / fitted_params[label])
            abs_errors_dev.append(fitted_params[label] - true_params[label])
            true_vert_dev.append(true_params[label])  
        else:
            val_norm_errors_drift_vel.append((fitted_params[label] - true_params[label]) / fitted_params[label])
            abs_errors_vel.append(fitted_params[label] - true_params[label])
            true_drift_vel.append(true_params[label])


print ""


# Plot absolute differences between fitted, true parameters

#
## Plot Histograms of true parameter values for detector, normalised differences between fitted, true values. Also plot expected pdfs for these values. 
#

# Plot true parameter values
bin_contents, bin_edges, _ = plt.hist(true_drift_vel, 20)
bin_width = (bin_edges[-1] - bin_edges[0]) / (len(bin_edges) - 1)
x_gaus = np.linspace(bin_edges[0], bin_edges[-1], 500)
plt.plot(x_gaus, sum(bin_contents) * bin_width * mlab.normpdf(x_gaus, 0, 0.02), 'r')
plt.title("True Values of Drift Velocity Deviation for Detector")
plt.ylabel("Number of Detector Planes")
plt.xlabel("Drift Velocity Deviation")
plt.show()
plt.clf()

bin_contents, bin_edges, _ = plt.hist(true_vert_dev, 20)
bin_width = (bin_edges[-1] - bin_edges[0]) / (len(bin_edges) - 1)
x_gaus = np.linspace(bin_edges[0], bin_edges[-1], 500)
plt.plot(x_gaus, sum(bin_contents) * bin_width * mlab.normpdf(x_gaus, 0, 0.1), 'r')
plt.title("True Values of Plane Displacement Deviation for Detector")
plt.ylabel("Plane Count")
plt.xlabel("Plane Displacement Deviation")
plt.show()
plt.clf()

if (len(fitted_param_errors) > 0):
    # Plot differences in true, fitted parameters normalised by fitted parameter uncertainty
    bin_contents, bin_edges, _ = plt.hist(unc_norm_errors_drift_vel, 20)
    bin_width = (bin_edges[-1] - bin_edges[0]) / (len(bin_edges) - 1)
    x_gaus = np.linspace(bin_edges[0], bin_edges[-1], 500)
    plt.plot(x_gaus, sum(bin_contents) * bin_width * mlab.normpdf(x_gaus, 0, 1), 'r')
    plt.title("Differences Between Fitted, True Values of DVD for Detector Planes, \n Normalised by Fitted Parameter Uncertainty")
    plt.ylabel("Number of Detector Planes")
    plt.xlabel("Normalised Parameter Difference")
    plt.show()
    plt.clf()

    bin_contents, bin_edges, _ = plt.hist(unc_norm_errors_vert_dev, 20)
    bin_width = (bin_edges[-1] - bin_edges[0]) / (len(bin_edges) - 1)
    x_gaus = np.linspace(bin_edges[0], bin_edges[-1], 500)
    plt.plot(x_gaus, sum(bin_contents) * bin_width * mlab.normpdf(x_gaus, 0, 1), 'r')
    plt.title("Differences Between Fitted, True Values of PDD for Detector Planes, \n Normalised by Fitted Parameter Uncertainty")
    plt.xlabel("Normalised Parameter Difference")
    plt.ylabel("Number of Detector Planes")
    plt.show()
    plt.clf()


# Plot differences in true, fitted parameters normalised by fitted parameter value
plt.hist(val_norm_errors_drift_vel, 20)
plt.title("Differences Between Fitted, True Values of DVD for Detector Planes, \n Normalised by Fitted Parameter Value")
plt.ylabel("Number of Detector Planes")
plt.xlabel("Normalised Parameter Difference")
plt.show()
plt.clf()

plt.hist(val_norm_errors_vert_dev, 20)
plt.title("Differences Between Fitted, True Values of PDD for Detector Planes, \n Normalised by Fitted Parameter Value")
plt.xlabel("Normalised Parameter Difference")
plt.ylabel("Number of Detector Planes")
plt.show()
plt.clf()


# Plot absolute differences between fitted, true parameters
plt.hist(abs_errors_vel, 20)
plt.title("Differences Between Fitted, True Values of DVD for Detector Planes")
plt.ylabel("Number of Detector Planes")
plt.xlabel("Fitted, True Parameter Difference")
plt.show()
plt.clf()

plt.hist(abs_errors_dev, 20)
plt.title("Differences Between Fitted, True Values of PDD for Detector Planes")
plt.ylabel("Number of Detector Planes")
plt.xlabel("Fitted, True Parameter Difference")
plt.show()
plt.clf()
