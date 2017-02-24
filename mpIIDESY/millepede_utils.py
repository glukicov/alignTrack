#
# Module to assist in saving parameters to a root file.
#
# John Smeaton 25/01/2017
#

import ROOT
from rootpy.io import root_open
from rootpy.tree import Tree, TreeModel
from rootpy.tree import IntCol, DoubleCol
import os
import shutil

# Class to create tree of parameters
class ParamTree(TreeModel):
    label = IntCol()
    fitType = IntCol()
    paramValue = DoubleCol()
    paramError = DoubleCol()


# Class for saving parameters from mille and pede to root files.
class MillepedeRootSaving: 

    # Takes arguments of filename for fitted parameters, and for steering file.
    def __init__(self, fitted_params_filename="millepede.res", steering_filename="mp2test1str_c.txt", true_params_filename="mp2test1_true_params_c.txt", root_filename="mptest1_parameters_c.root"):
        
        # Set variables from arguments
        self.fitted_params_filename = fitted_params_filename
        self.true_params_filename = true_params_filename
        self.steering_filename = steering_filename
        self.root_filename = root_filename

        # Dictionary of fit types
        self.fit_type_dict = {1: "inversion", 2: "diagonalization", 3: "fullMINRES", 4: "sparseMINRES"}


    # Function to create an empty tree of parameter values
    def create_tree(self):
        f = root_open(self.root_filename, "recreate")
        t = Tree("paramTree", model=ParamTree)
        t.write("", ROOT.TObject.kWriteDelete)
        f.close()


    # Function to save parameters from a pede output file to root
    def save_to_root(self, fit_type):

        # Open root file, and create tree
        f = root_open(self.root_filename, "update")  
        t = f.paramTree

        # Cannot write to tree buffer without this, for some reason.
        for entries in t:
            continue

        # Check if saving true or fitted parameter values
        if (fit_type == 0):
            params_f = open(self.true_params_filename, 'r')
        else:
            params_f = open(self.fitted_params_filename, 'r')

        for line in params_f.readlines():
        
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


    # Function to set whether steering file uses a given input file. Takes arguments of boolean for whether to read file, and what filename to toggle.
    def set_input_file_reading(self, to_read, filename_to_toggle):
        
        # New steering file to write lines to
        new_steering_f = open(self.steering_filename + ".tmp", "w")
        
        # Iterate across lines of original steering file
        with open(self.steering_filename, 'r') as steering_f:
            for line in steering_f.readlines():

                # Split line into items, delimiting by spaces
                items = line.split()

                # Skip line if empty
                if (len(items)) == 0:
                    new_steering_f.write(line)
                    continue

                # Check if this line is the filename to toggle. Comment out line if don't want to read this input file.
                if ((items[0] == filename_to_toggle) and not (to_read)):
                    items[0] = "*" + filename_to_toggle
                    newline = ' '.join(items) + "\n"
                    new_steering_f.write(newline)
                    
                # Check if this line is filename commented out. Uncomment if want to read this input file.
                elif ((items[0] == "*" + filename_to_toggle) and to_read):
                    items[0] = filename_to_toggle
                    newline = ' '.join(items) + "\n"
                    new_steering_f.write(newline)

                # Otherwise, write line to new steering file
                else:
                    new_steering_f.write(line)

        # Move new steering file to old steering filename.
        shutil.move(self.steering_filename + ".tmp", self.steering_filename)


    # Function to modify the pede steering file, so it carries out the specified fit type.
    def modify_steering_file_fit(self, fit_type):
        
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


    # Function to remove a given constraint from the constraint file.
    def remove_constraint(self, constraint_filename, constraint_number):

        # New constraint file
        new_constraint_file = open(constraint_filename + ".tmp", "w")

        # Number of constraints previously seen in file
        constraints_seen = 0

        with open(constraint_filename, "r") as constraint_file:
            for line in constraint_file:

                items = line.split()

                # Check line long enough to process. If not, write old line to new file
                if (len(items) < 0):
                    new_constraint_file.write(line)
                    continue

                # If line corresponds to new constraint, increment number of constraints seen
                if (items[0] == "Constraint"):
                    constraints_seen += 1

                # If number of constraints seen corresponds to number of constraints previously seen, don't copy this line
                if not (constraints_seen == constraint_number):
                    new_constraint_file.write(line)

        # Move new constraint file to old constraint filename.
        shutil.move(constraint_filename + ".tmp", constraint_filename)
