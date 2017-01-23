import matplotlib.pyplot as plt 
import numpy as np
import sys
import getopt


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
errors_drift_vel = []
errors_vert_dev = []

# Iterate across all keys in dictionary
for label in sorted(true_params.iterkeys()):

    # Print blank line between PDD, DVD outputs 
    if label == 501: print ""

    # Print parameter label, true and fitted values, fitted uncertainty, and normalised difference between fitted, true values
    print label, "\t", true_params[label], "\t", fitted_params[label], "\t", fitted_param_errors[label], "\t", (fitted_params[label] - true_params[label]) / fitted_param_errors[label]

    # Add normalised difference, true parameter values to appropriate lists, checking label value to see if entry corrsponds to displacement of velocity deviation
    if label < 500:
        errors_vert_dev.append((fitted_params[label] - true_params[label]) / fitted_param_errors[label])
        true_vert_dev.append(true_params[label])  
    else:
        errors_drift_vel.append((fitted_params[label] - true_params[label]) / fitted_param_errors[label])
        true_drift_vel.append(true_params[label])


#
## Plot Histograms of true parameter values for detector, differences between fitted, true values.
#

plt.hist(true_drift_vel, 20)
plt.title("True Values of Drift Velocity Deviation for Detector")
plt.ylabel("Number of Detector Planes")
plt.xlabel("Drift Velocity Deviation")
plt.show()
plt.clf()

plt.hist(true_vert_dev, 20)
plt.title("True Values of Plane Displacement Deviation for Detector")
plt.ylabel("Plane Count")
plt.xlabel("Plane Displacement Deviation")
plt.show()
plt.clf()


plt.hist(errors_drift_vel, 20)
plt.title("Differences Between Fitted, True Values of DVD for Detector Planes, \n Normalised by Fitted Parameter Uncertainty")
plt.ylabel("Number of Detector Planes")
plt.show()
plt.clf()

plt.hist(errors_vert_dev, 20)
plt.title("Differences Between Fitted, True Values of PDD for Detector Planes, \n Normalised by Fitted Parameter Uncertainty")
plt.ylabel("Number of Detector Planes")
plt.show()
plt.clf()
