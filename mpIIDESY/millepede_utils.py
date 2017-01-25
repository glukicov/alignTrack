#
# Module to assist in saving parameters to a root file.
#
# John Smeaton 25/01/2017
#

import ROOT
from rootpy.io import root_open
import os
import shutil


# Class for saving parameters from mille and pede to root files.
class MillepedeRootSaving: 

    # Takes arguments of filename for fitted parameters, and for steering file.
    def __init__(self, fitted_params_filename="millepede.res", steering_filename="mp2test1str_c.txt"):
        
        # Set variables from arguments
        self.fitted_params_filename = fitted_params_filename
        self.steering_filename = steering_filename
        self.root_filename = "mptest1_parameters.root"

        # Dictionary of fit types
        self.fit_type_dict = {1: "inversion", 2: "diagonalization", 3: "fullMINRES", 4: "sparseMINRES"}


    # Function to save parameters from a pede output file to root
    def save_to_root(self, fit_type):

        # Open root file, and create tree
        f = root_open(self.root_filename, "update")  
        t = f.paramTree

        # Cannot write to tree buffer without this, for some reason.
        for entries in t:
            continue

        # Open text file of fitted parameters, iterating across all lines
        with open(self.fitted_params_filename, 'r') as fitted_f:
            for line in fitted_f.readlines():
        
                # Split line into items 
                items = line.split()        

                # Test if first entry in file is a label (can cast to int). If not, continue.
                try:
                    int(items[0])
                except:
                    continue

                # Save type of fit, and label to tree
                t.fitType = fit_type;
                t.label = int(items[0])
        
                # Save fitted parameter value to tree. If it doesn't exist, continue.
                try:
                    t.paramValue = float(items[1])
                except Exception:
                    continue

                # Save parameter error, if it exists. If not, set to zero.
                try:
                    t.paramError = float(items[4])
                except Exception:
                    t.paramError = 0.0

                # Fill tree
                t.fill();
                    
        # Write to tree, and close file
        t.write("", ROOT.TObject.kWriteDelete)
        f.close()


    # Function to modify the pede steering file, so it carries out the specified fit type.
    def modify_steering_file(self, fit_type):
        
        # New steering file to write lines to
        new_steering_f = open(self.steering_filename + ".tmp", "w")
        
        # Iterate across lines of original steering file
        with open(self.steering_filename, 'r') as steering_f:
            for line in steering_f.readlines():

                # Split line into items, delimiting by spaces
                items = line.split()
            
                # Check line long enough to process. If not, write old line to new file
                if (len(items) < 2):
                    new_steering_f.write(line)
                    continue
            
                # If line is for method, and this method is not the desired one, set steering file so this line is ignored.
                if (items[0] == "method") and not (items[1] == self.fit_type_dict[fit_type]):
                    items[0] = "*method"
                    newline = ' '.join(items) + "\n"
                    new_steering_f.write(newline)

                # If line is for desired method, and line is set to be ignored, set so line is no longer ignored.
                elif (items[0] == "*method") and (items[1] == self.fit_type_dict[fit_type]):
                    items[0] = "method"
                    newline = ' '.join(items) + "\n"
                    new_steering_f.write(newline)
                
                # If line is for anything else, replicate in new steering file.
                else:
                    new_steering_f.write(line)
                    
        # Move new steering file to old steering filename.
        shutil.move(self.steering_filename + ".tmp", self.steering_filename)



