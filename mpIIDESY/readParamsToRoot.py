#
# Script to generate data using the C++ port of MpTest1, then carry out fits 
# with pede, using various methods, and saving true and fitted parameters to 
# a Root TTree.
#
# John Smeaton 25/01/2017
#

import millepede_utils
import os

fitted_params_filename="millepede.res"
steering_filename="mp2test1str_c.txt"

# Set up instance of utilities class
m_utils = millepede_utils.MillepedeRootSaving()

# Generate mille data.
os.system("./MpTest1")

# Iterate across types of fit.
for i in m_utils.fit_type_dict.iterkeys():

    # Change steering file to carry out this type of fit, then run pede and save results to root.
    m_utils.modify_steering_file(i)
    os.system("./pede mp2test1str_c.txt")
    m_utils.save_to_root(i)
