#
# Script to generate data using the original Fortran version of MpTest1, then carry out fits with pede, using various methods, and saving true and fitted 
# parameters to a Root TTree.
#
# John Smeaton 07/02/2017
#

import millepede_utils
import os

# Names of filenames for fitted parameter results, and steering file
fitted_params_filename = "millepede.res"
true_params_filename = "mp2test1_true_params_fortran.txt"
steering_filename = "mp2str.txt"
root_filename = "mptest1_parameters_fortran.root"

# Set up instance of utilities class
m_utils = millepede_utils.MillepedeRootSaving(steering_filename, true_params_filename, fitted_params_filename, root_filename)

# Remove old versions of steering file, root output
try:
    os.system("rm " + fitted_params_filename)
    os.system("rm " + true_params_filename)
    os.system("rm " + root_filename)
    os.system("rm " + steering_filename)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise

# Create a new root file with an empty tree for parameter values
m_utils.create_tree() 

# Run MpTest1 Fortran version, using pede test flag, generating mille data
os.system("./pede -t")

# Save true parameter values to root tree.
m_utils.save_to_root(0)

# Iterate across types of fit.
for i in m_utils.fit_type_dict.iterkeys():

    # Change steering file to carry out this type of fit, then run pede and save results to root.
    m_utils.modify_steering_file_fit(i)
    os.system("./pede " + steering_filename)
    m_utils.save_to_root(i)
