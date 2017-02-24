""" Script to perform eigenvalue analysis of solution, in effort to remove 
weak modes. Opens file of eigenvectors output by "Diagonalisation" method of
Millepede solution, reading values for analysis.

Author: John Smeaton
Date: 21/02/2017"""


import sys
import getopt
import os.path
import numpy as np
import matplotlib.pyplot as plt


class Eigenvector():
    """Class to represent a single eigenvector in the solution"""

    def __init__(self, eigenvec_index, eigenval):
        # Set values for eigenvalue, index of this vector, and initialise
        # dictionary of eigenvector values for each parameter.
        self.eigenval = eigenval
        self.eigenvec_index = eigenvec_index
        self.param_vals = dict()

    def __str__(self):
        multiline_string = ("Eigenvector " + str(self.eigenvec_index) + 
                            " with eigenvalue " + str(self.eigenval) + "\n")
        for param_index in sorted(self.param_vals):
            multiline_string += ("Param. " + str(param_index) + 
                                 " " + str(self.param_vals[param_index]) + 
                                 "\n")
        return multiline_string


# Get system arguments, define string showing help message
argv = sys.argv[1:]
helpstring = "eigenSpectrumAnalysis.py -i <input_file> -o <output_dir>"

# Define input filename, output directory strings
input_filename = ""
output_dir = ""

# Get options, arguments
try:
    opts, args = getopt.getopt(argv, "h:i:o:", ["help", "input_file", "output_dir"])
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
    elif opt in ("-o", "--output_dir"):
        output_dir = os.path.abspath(arg)

# Creat output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# Eigenvectors found in file, and variable for eigenvector being currently 
# read
found_eigenvecs = []
current_eigenvec = None

# Open input file, and iterate over lines
with open(input_filename) as eigenval_output:
    for line in eigenval_output.readlines():

        # Split lines into items
        items = line.split()

        # Skip this line if it's empty
        if (len(items) == 0):
            continue

        # Start new eigenvector if this line marked as such, and append 
        # previous eigenvector to list if appropriate
        if (items[0] == "Eigenvector"):
            if current_eigenvec is not None:
                found_eigenvecs.append(current_eigenvec)
            current_eigenvec = Eigenvector(int(items[1]), float(items[4]))
            continue

        # Current index of parameter, and its value in eigenvector
        current_param_index = 0
        current_param_val = 0

        # Check start of line is an integer (should always be case, but just 
        # in case.)
        if not (float(items[0]).is_integer()):
            continue

        # Iterate over items in line
        for item in items:
            
            # Get current parameter index if this is an integer, or parameter
            # value if not, then insert current parameter value and index 
            # into dictionary.  
            if (float(item).is_integer() and float(item) > 0.0):
                current_param_index = int(item)
            else:
                current_param_val = float(item)
                current_eigenvec.param_vals[current_param_index] = current_param_val

    # After last line, append eigenvector to list.
    found_eigenvecs.append(current_eigenvec)


# Get indices for eigenvectors
eigenvec_indeces = [str(eigenvec.eigenvec_index) for eigenvec in found_eigenvecs]
eigenvals = [eigenvec.eigenval for eigenvec in found_eigenvecs]

# Where to show bars, and width
ind = np.arange(len(found_eigenvecs))
width = 1.0

# Plot bar chart of eigenvalues
fig, ax = plt.subplots()
ax.bar(ind, eigenvals, width)

# Set axis scaling and tick labels on x-axis
ax.set_yscale('log')
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(eigenvec_indeces)

# Set titles
ax.set_title("Eigenspectrum for System")
ax.set_xlabel("Eigenvector Index")
ax.set_ylabel("Eigenvalue")

plt.savefig(output_dir + "/eigenspectrum.pdf")
plt.clf()


solution_vals = dict()

for eigenvec in found_eigenvecs:

    label_plane_disp = []
    fit_plane_disp = []

    label_vel_dev = []
    fit_vel_dev = []

    for param_label in eigenvec.param_vals:

        if param_label < 500:
            label_plane_disp.append(param_label)
            fit_plane_disp.append(eigenvec.param_vals[param_label])
        else:
            label_vel_dev.append(param_label)
            fit_vel_dev.append(eigenvec.param_vals[param_label])


    # Plot scatter graph of plane displacements against parameter label
    plt.plot(label_plane_disp, fit_plane_disp, 'b.')
    plt.title("Fitted Plane Displacement Parameter Values \n Eigenvector " + str(eigenvec.eigenvec_index) + " with Eigenvalue " + str(eigenvec.eigenval))
    plt.ylabel("Fitted Plane Displacement")
    plt.xlabel("Parameter Label")
    plt.savefig(output_dir + "/pdd_eigenvec_" + str(eigenvec.eigenvec_index) + ".pdf")
    plt.clf()

    # Plot scatter graph of plane displacements against parameter label
    plt.plot(label_vel_dev, fit_vel_dev, 'b.')
    plt.title("Fitted Velocity Deviation Parameter Values \n Eigenvector " + str(eigenvec.eigenvec_index) + " with Eigenvalue " + str(eigenvec.eigenval))
    plt.ylabel("Fitted Velocity Deviation")
    plt.xlabel("Parameter Label")
    plt.savefig(output_dir + "/dvd_eigenvec_" + str(eigenvec.eigenvec_index) + ".pdf")
    plt.clf()
