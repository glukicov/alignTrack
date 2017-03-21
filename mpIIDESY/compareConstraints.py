#
# Script which creates data using MpTest1, then carries out pede fits, varying 
# the constraints used. The differences between the fitted and true plane 
# displacements are then plotted.
#
# John Smeaton 24/02/2017
#

import millepede_utils
import os
import numpy as np
import matplotlib.pyplot as plt


# Function to read plane deviations from a specified filename
def read_plane_deviations(file_to_read):

    # List of plane deviations read from file
    plane_deviations = []

    # Open file, and iterate over lines
    with open(file_to_read, "r") as params_f:
        for line in params_f.readlines():
        
            # Split line into items 
            items = line.split()        

            # Test if first entry in file is a label (can cast to int). If not, continue.
            try:
                int(items[0])
            except:
                continue

            # Test if label for parameter represents a plane deviation, rather than a velocity deviation
            if int(items[0]) < 500:
                try:
                    plane_deviations.append(float(items[1]))
                except Exception:
                    continue
        
    # Cast list to a numpy array, and return
    return np.array(plane_deviations)


# Function to read errors on fitted plane deviations from a specified filename
def read_plane_errors(file_to_read):

    # List of plane errors read from file
    plane_errors = []

    # Open file, and iterate over lines
    with open(file_to_read, "r") as params_f:
        for line in params_f.readlines():
        
            # Split line into items 
            items = line.split()        

            # Test if first entry in file is a label (can cast to int). If not, continue.
            try:
                int(items[0])
            except:
                continue

            # Test if label for parameter represents a plane deviation, rather than a velocity deviation
            if int(items[0]) < 500:
                try:
                    plane_errors.append(float(items[4]))
                except Exception:
                    continue
        
    # Cast list to a numpy array, and return
    return np.array(plane_errors)



# Names of filenames for fitted parameter results, and steering file
fitted_params_filename = "millepede.res"
true_params_filename="mp2test1_true_params_c.txt"
steering_filename = "mp2test1str_c.txt"
constraint_filename = "mp2test1con_c.txt"
root_filename = "mptest1_parameters_c.root"

# Set up instance of utilities class
m_utils = millepede_utils.MillepedeRootSaving(root_filename=root_filename, steering_filename=steering_filename, true_params_filename=true_params_filename, fitted_params_filename=fitted_params_filename)

# Remove old versions of steering file, root output
try:
    os.system("rm mp2str.txt")
    os.system("rm mptest1_parameters_c.root")
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise

# Generate mille data.
os.system("./MpTest1 uniform_ran.txt gaussian_ran.txt")

#
# All constraints
#
# Read true plane deviations
plane_deviations_true = read_plane_deviations(true_params_filename)

# Carry out fit, then read fitted plane deviations
os.system("./pede " + steering_filename)
plane_deviations_all_cons = read_plane_deviations(fitted_params_filename)
plane_errors_all_cons = read_plane_errors(fitted_params_filename)


#
# No shift constraint
#
# Remove plane shift constraint from file
m_utils.remove_constraint(constraint_filename, 1)

# Carry out fit, then read fitted plane deviations
os.system("./pede " + steering_filename)
plane_deviations_no_shift_cons = read_plane_deviations(fitted_params_filename)
plane_errors_no_shift_cons = read_plane_errors(fitted_params_filename)


#
# No shear constraint
#
# Recreate constraint file, then remove plane shear constraint
os.system("./MpTest1 uniform_ran.txt gaussian_ran.txt")
m_utils.remove_constraint(constraint_filename, 2)

# Carry out fit, then read fitted plane deviations
os.system("./pede " + steering_filename)
plane_deviations_no_shear_cons = read_plane_deviations(fitted_params_filename)
plane_errors_no_shear_cons = read_plane_errors(fitted_params_filename)


#
# No constraints
#
# Remove constraint file from steering file, carry out fit, and read fitted plane deviations
m_utils.set_input_file_reading(False, constraint_filename)
os.system("./pede " + steering_filename)
plane_deviations_no_cons = read_plane_deviations(fitted_params_filename)
plane_errors_no_cons = read_plane_errors(fitted_params_filename)

# Generate array of plane displacement labels
labels = np.arange(12, 212, 2)


lower_tolerance = np.array([-0.015 for i in labels])
upper_tolerance = np.array([0.015 for i in labels])

# Plot parameter differences
plt.fill_between(labels, lower_tolerance, upper_tolerance, color="yellow", alpha=0.5, label="1$\sigma$ Tolerance")
plt.errorbar(labels, plane_deviations_all_cons - plane_deviations_true, color='b', fmt='.', label="All Constraints", xerr=0, yerr=plane_errors_all_cons)
plt.errorbar(labels, plane_deviations_no_shift_cons - plane_deviations_true, color='g', fmt='.', label="No Shift Constraint", xerr=0, yerr=plane_errors_no_shift_cons)
plt.errorbar(labels, plane_deviations_no_shear_cons - plane_deviations_true, color='m', fmt='.', label="No Shear Constraint", xerr=0, yerr=plane_errors_no_shear_cons)
plt.errorbar(labels, plane_deviations_no_cons - plane_deviations_true, color='r', fmt='.', label="No Constraints", xerr=0, yerr=plane_errors_no_cons)

plt.xlim([12,210])

# Label plot
plt.title("Unnormalised Differences Between Fitted, True \n Plane Displacements for Different Constraints")
plt.ylabel("Difference Between Fitted, True Plane Displacement")
plt.xlabel("Parameter Label")
plt.legend(loc='best')

plt.show()

