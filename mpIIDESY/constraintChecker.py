"""Script to check constraints on parameter values are obeyed, using pede results file and constraints file. Currently checks overall detector dispacement and shear. Outputs sum of constraint weightings multiplied by parameter results. If these are not equal to zero (or far off), the constraint is not obeyed.

Author: John Smeaton
Date: 27/02/2017"""

import sys
import getopt
import os.path
import numpy as np
import matplotlib.pyplot as plt


# Get system arguments, define string showing help message
argv = sys.argv[1:]
helpstring = "constraintChecker.py -r <results_file> -c <constraints_file>"

# Define input filename, output directory strings
result_filename = ""
constraint_filename = ""

# Get options, arguments
try:
    opts, args = getopt.getopt(argv, "h:r:c:", ["help", "results_file", "constraints_file"])
except getopt.GetoptError:
    print helpstring
    sys.exit(2)

# Set input filename, if given in console argument
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print helpstring
        sys.exit(2)
    elif opt in ("-r", "--results_file"):
        result_filename = os.path.abspath(arg)
    elif opt in ("-c", "--constraints_file"):
        constraint_filename = os.path.abspath(arg)


# Containers for parameter results, and constraints on parameters
param_results = dict()
param_constraints = []

# Open results file, and iterate over lines
with open(result_filename) as result_file:
    for line in result_file.readlines():

        # Split line into items
        items = line.split()

        # Skip this line if it's empty
        if (len(items) == 0):
            continue

        # Add result to dictionary
        try:
            int(items[0])
            param_results[int(items[0])] = float(items[1])
        except Exception:
            pass


# Open constraints file, and iterate over lines.
with open(constraint_filename) as constraint_file:
    for line in constraint_file.readlines():

        # Split line into items
        items = line.split()

        # Skip this line if empty
        if (len(items) == 0):
            continue

        # Add new dictionary of constraint weightings, if this line shows a new constraint
        if(items[0] == "Constraint"):
            param_constraints.append(dict())

        # Add constraint weighting to current dictionary of constraints
        try:
            int(items[0])
            param_constraints[-1][int(items[0])] = float(items[1])
        except Exception:
            pass


# Iterate over each constraint
for constraint in param_constraints:

    constraint_sum = 0 # Sum of parameter values multiplied by constraint weightings

    # Add to sum
    for label in constraint:
        constraint_sum += constraint[label] * param_results[label]

    # Print sum
    print "Constraint Sum:", constraint_sum 
