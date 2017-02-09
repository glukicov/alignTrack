import toytraceback as ttb
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys
import math
import scipy.stats as stats
import os.path
import getopt
import ROOT
import psutil
import os
import subprocess

# Script to produce scatter plots of module misalignment against chi-squared, to examine correlations

# Binds script to single CPU (Prevents problems on batch farm)
psu=psutil.Process()
mypid=os.getpid()
pscommand=subprocess.Popen("/bin/ps -o psr "+str(mypid),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
psout=pscommand.stdout.read()
pscommand.wait()
psu.cpu_affinity([int(psout.split('\n')[1])])

# Root file input, and output directory for scatter plots 
input_file = ""
output_dir = ""

# Cut on chi-squared to examine plots above and below, and number of points to show in scatter plot
chi_sq_cut = 2
points_shown = 250

# Get arguments
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "h:i:o:", ["help", "input_file=", "output_dir"])
except getopt.GetoptError:
    print "tree_reader.py -i <input_file> -o <output_dir>"
    sys.exit(2)

# Set parameters according to arguments
for opt, arg in opts:

    if opt in ("-h", "--help"):
        print "tree_reader.py -i <input_file> -o <output_dir>"
        sys.exit()
    elif opt in ("-i", "--input_file"):
        input_file = arg
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    

print input_file

# Get input file
f = ROOT.TFile.Open(input_file, "read")

# Get tree of straw hit events
tree = f.Get("event_tree")

# Lists of module alignments, and dictionary of lists of reduced chi-squared values for different degreees of freedom
module_misalignments = dict()
chi_sq_vals = dict()

# Get number of fits in tree, by examining last entry in tree 
entry_count = tree.GetEntries()
tree.GetEntry(entry_count - 1)
fit_count = tree.fitNum + 1

track_count = 0 # For number of tracks of fits in file

print "Total Fits:", str(fit_count)

# Loop across all fits in file
for i in xrange(fit_count):
    
    if (i % 100) == 0:
        print i
        #memtrack.print_diff()

    # Get entries with correct fit number in tree
    fit_cut_tree = tree.CopyTree('fitNum==' + str(i))

    # Get track count for this fit
    fit_cut_tree.GetEntry(fit_cut_tree.GetEntries() - 1)
    track_count = fit_cut_tree.eventNum + 1

    # Lists of hit distances for fit
    true_hit_distances = []
    fitted_hit_distances = []
 
    # Get number of entries for this fit, and loop over them
    fit_cut_tree_entry_count = fit_cut_tree.GetEntries() 
            
    for j in xrange(fit_cut_tree_entry_count):
        
        fit_cut_tree.GetEntry(j)

        # Append hit distances for this entry
        true_hit_distances.append(fit_cut_tree.trueHitDistance)
        fitted_hit_distances.append(fit_cut_tree.fittedHitDistance)

    # Create array of hit distance uncertainties, taken to be equal to hit resolution of detector    
    sigmas = np.array([ttb.hit_resolution for n in xrange(fit_cut_tree_entry_count)])

    # Calculate reduced value of chi-squared from hit distances
    chi_sq_red = np.sum(((np.array(true_hit_distances) - np.array(fitted_hit_distances)) / sigmas)**2) / (fit_cut_tree_entry_count - ((2 * track_count) + 1)) 

    fit_cut_tree.GetEntry(0)
    
    # Add misalignment to dictionary
    if fit_cut_tree_entry_count - ((2 * track_count) + 1) in module_misalignments:
        module_misalignments[fit_cut_tree_entry_count - ((2 * track_count) + 1)].append(fit_cut_tree.fittedAlignment - fit_cut_tree.trueAlignment)
    else:
        module_misalignments[fit_cut_tree_entry_count - ((2 * track_count) + 1)] = [fit_cut_tree.fittedAlignment - fit_cut_tree.trueAlignment]
            

    # Add chi-squared value to dictionary
    if fit_cut_tree_entry_count - ((2 * track_count) + 1) in chi_sq_vals:
        chi_sq_vals[fit_cut_tree_entry_count - ((2 * track_count) + 1)].append(chi_sq_red)
    else:
        chi_sq_vals[fit_cut_tree_entry_count - ((2 * track_count) + 1)] = [chi_sq_red]

    fit_cut_tree.IsA().Destructor(fit_cut_tree) # To avoid stupid bloody memory leak

# Close input file
f.Close()

# Arrays to contain values for all degrees of freedom
chi_sq_vals_alldof = []
module_misalignments_alldof = []

# Generate string describing whether hit distance limit is in place
straw_status = ""
if len(chi_sq_vals) > 1:
    straw_status = "missed_straws"
else:
    straw_status = "all_straws"

# Loop over degrees of freedom of recorded chi-squared values
for dof in chi_sq_vals:

    # Add entries to array with all degrees of freedom
    chi_sq_vals_alldof.extend(chi_sq_vals[dof])
    module_misalignments_alldof.extend(module_misalignments[dof])

    # Get sample of values from array. If less points than desired, set sample to array of all points
    chi_sq_vals_sample = []
    module_misalignments_sample = []
    if len(module_misalignments[dof]) >= points_shown:
        module_misalignments_sample = module_misalignments[dof][:points_shown]
        chi_sq_vals_sample = chi_sq_vals[dof][:points_shown]
    else:
        module_misalignments_sample = module_misalignments[dof]
        chi_sq_vals_sample = chi_sq_vals[dof]
        
    # Array of values above and below chi-squared cuts
    chi_sq_vals_abovecut = []
    chi_sq_vals_belowcut = []
    module_misalignments_abovecut = []
    module_misalignments_belowcut = []
    
    # Populate arrays with values above, below cut
    for i in xrange(len(chi_sq_vals[dof])):
        if chi_sq_vals[dof][i] > chi_sq_cut:
            chi_sq_vals_abovecut.append(chi_sq_vals[dof][i])
            module_misalignments_abovecut.append(module_misalignments[dof][i])
        else:
            chi_sq_vals_belowcut.append(chi_sq_vals[dof][i])
            module_misalignments_belowcut.append(module_misalignments[dof][i])
            
    # Generate matrices of covariance with and without cut, and for misalignments and modulus of misalignments 
    cov_mat = np.cov([module_misalignments[dof], chi_sq_vals[dof]])
    cov_mat_absolute = np.cov([np.absolute(module_misalignments[dof]), chi_sq_vals[dof]])
    cov_mat_abovecut = np.cov([module_misalignments_abovecut, chi_sq_vals_abovecut])
    cov_mat_absolute_abovecut = np.cov([np.absolute(module_misalignments_abovecut), chi_sq_vals_abovecut])
    cov_mat_belowcut = np.cov([module_misalignments_belowcut, chi_sq_vals_belowcut])
    cov_mat_absolute_belowcut = np.cov([np.absolute(module_misalignments_belowcut), chi_sq_vals_belowcut])


    print ""
    print ""
    print str(dof), " Degrees of Freedom"
    print ""
    print "No Cut"
    print ""
    print "Misalignment and Chi-Squared:"
    print "Covariance Matrix"
    print cov_mat
    print ""
    print "Correlation: ", str(cov_mat[0][1] / math.sqrt(cov_mat[0][0] * cov_mat[1][1]))
    print ""
    print ""

    # Plot scatter, and label
    plt.plot(module_misalignments_sample, chi_sq_vals_sample, 'bo')
    plt.title("Plot of $\chi^2_{red}$ Against Module Misalignment, \n for Each Fit, With " + str(dof) + " D.o.f.")
    plt.ylabel("$\chi^2_{red}$")
    plt.xlabel("Module Misalignment / mm")

    # Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
    number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals[dof]))))
    shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len(chi_sq_vals_sample))))
    correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat[0][1] / math.sqrt(cov_mat[0][0] * cov_mat[1][1]))))

    # Get handles and labels from patches, and create legend
    handles, labels = plt.subplot(111).get_legend_handles_labels()
    handles.append(number_patch)
    handles.append(shown_patch)
    handles.append(correlation_patch)
    labels.append(number_patch.get_label())
    labels.append(shown_patch.get_label())
    labels.append(correlation_patch.get_label())
    plt.legend(handles, labels, loc='best', frameon=False)

    # Save figure, and clear
    plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_" +str(dof) + "dof_nocut.pdf")
    plt.clf()


    print "Absolute Misalignment and Chi-Squared"
    print "Covariance Matrix (Absolute Misalignment)"
    print cov_mat_absolute
    print ""
    print "Correlation: ", str(cov_mat_absolute[0][1] / math.sqrt(cov_mat_absolute[0][0] * cov_mat_absolute[1][1]))

    # Plot scatter, and label
    plt.plot(np.absolute(module_misalignments_sample), chi_sq_vals_sample, 'bo')
    plt.title("Plot of $\chi^2_{red}$ Against Modulus of Module Misalignment, \n for Each Fit, With " + str(dof) + " D.o.f.")
    plt.ylabel("$\chi^2_{red}$")
    plt.xlabel("Absolute Module Misalignment / mm")

    # Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
    number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals[dof]))))
    shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len(chi_sq_vals_sample))))
    correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_absolute[0][1] / math.sqrt(cov_mat_absolute[0][0] * cov_mat_absolute[1][1]))))

    # Get handles and labels from patches, and create legend
    handles, labels = plt.subplot(111).get_legend_handles_labels()
    handles.append(number_patch)
    handles.append(shown_patch)
    handles.append(correlation_patch)
    labels.append(number_patch.get_label())
    labels.append(shown_patch.get_label())
    labels.append(correlation_patch.get_label())
    plt.legend(handles, labels, loc='best', frameon=False)

    # Save figure, and clear
    plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_" +str(dof) + "dof_nocut_absolute.pdf")
    plt.clf()


    print ""
    print ""
    print "Above Cut"
    print ""
    print "Misalignment and Chi-Squared:"
    print "Covariance Matrix"
    print cov_mat_abovecut
    print ""
    print "Correlation: ", str(cov_mat_abovecut[0][1] / math.sqrt(cov_mat_abovecut[0][0] * cov_mat_abovecut[1][1]))
    print ""
    print ""

    # Plot scatter, and label
    plt.plot(module_misalignments_sample, chi_sq_vals_sample, 'bo')
    plt.ylim(ymin=chi_sq_cut)
    plt.title("Plot of $\chi^2_{red}$ Against Module Misalignment, \n for Fits With $\chi^2_{red} > 2$, With " + str(dof) + " D.o.f.")
    plt.ylabel("$\chi^2_{red}$")
    plt.xlabel("Module Misalignment / mm")

    # Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
    number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_abovecut))))
    shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len([i for i in chi_sq_vals_sample if i > chi_sq_cut ]))))
    correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_abovecut[0][1] / math.sqrt(cov_mat_abovecut[0][0] * cov_mat_abovecut[1][1]))))

    # Get handles and labels from patches, and create legend
    handles, labels = plt.subplot(111).get_legend_handles_labels()
    handles.append(number_patch)
    handles.append(shown_patch)
    handles.append(correlation_patch)
    labels.append(number_patch.get_label())
    labels.append(shown_patch.get_label())
    labels.append(correlation_patch.get_label())
    plt.legend(handles, labels, loc='best', frameon=False)

    # Save figure, and clear
    plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_" +str(dof) + "dof_abovecut.pdf")
    plt.clf()


    print "Absolute Misalignment and Chi-Squared"
    print "Covariance Matrix (Absolute Misalignment)"
    print cov_mat_absolute_abovecut
    print ""
    print "Correlation: ", str(cov_mat_absolute_abovecut[0][1] / math.sqrt(cov_mat_absolute_abovecut[0][0] * cov_mat_absolute_abovecut[1][1]))

    # Plot scatter, and label
    plt.plot(np.absolute(module_misalignments_sample), chi_sq_vals_sample, 'bo')
    plt.title("Plot of $\chi^2_{red}$ Against Modulus of Module Misalignment, \n for Fits With $\chi^2_{red} > 2$, With " + str(dof) + " D.o.f.")
    plt.ylabel("$\chi^2_{red}$")
    plt.xlabel("Absolute Module Misalignment / mm")
    plt.ylim(ymin=chi_sq_cut)

    # Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
    number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_abovecut))))
    shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len([i for i in chi_sq_vals_sample if i > chi_sq_cut ]))))
    correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_absolute_abovecut[0][1] / math.sqrt(cov_mat_absolute_abovecut[0][0] * cov_mat_absolute_abovecut[1][1]))))

    # Get handles and labels from patches, and create legend
    handles, labels = plt.subplot(111).get_legend_handles_labels()
    handles.append(number_patch)
    handles.append(shown_patch)
    handles.append(correlation_patch)
    labels.append(number_patch.get_label())
    labels.append(shown_patch.get_label())
    labels.append(correlation_patch.get_label())
    plt.legend(handles, labels, loc='best', frameon=False)

    # Save figure, and clear
    plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_" +str(dof) + "dof_abovecut_absolute.pdf")
    plt.clf()


    print ""
    print ""
    print "Below Cut"
    print ""
    print "Misalignment and Chi-Squared:"
    print "Covariance Matrix"
    print cov_mat_belowcut
    print ""
    print "Correlation: ", str(cov_mat_belowcut[0][1] / math.sqrt(cov_mat_belowcut[0][0] * cov_mat_belowcut[1][1]))
    print ""
    print ""

    # Plot scatter, and label
    plt.plot(module_misalignments_sample, chi_sq_vals_sample, 'bo')
    plt.ylim(ymax=chi_sq_cut)
    plt.title("Plot of $\chi^2_{red}$ Against Module Misalignment, \n for Fits With $\chi^2_{red} < 2$, With " + str(dof) + " D.o.f.")
    plt.ylabel("$\chi^2_{red}$")
    plt.xlabel("Module Misalignment / mm")

    # Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
    number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_belowcut))))
    shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len([i for i in chi_sq_vals_sample if i < chi_sq_cut ]))))
    correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_belowcut[0][1] / math.sqrt(cov_mat_belowcut[0][0] * cov_mat_belowcut[1][1]))))

    # Get handles and labels from patches, and create legend
    handles, labels = plt.subplot(111).get_legend_handles_labels()
    handles.append(number_patch)
    handles.append(shown_patch)
    handles.append(correlation_patch)
    labels.append(number_patch.get_label())
    labels.append(shown_patch.get_label())
    labels.append(correlation_patch.get_label())
    plt.legend(handles, labels, loc='best', frameon=False)

    # Save figure, and clear
    plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_" +str(dof) + "dof_belowcut.pdf")
    plt.clf()


    print "Absolute Misalignment and Chi-Squared"
    print "Covariance Matrix (Absolute Misalignment)"
    print cov_mat_absolute_belowcut
    print ""
    print "Correlation: ", str(cov_mat_absolute_belowcut[0][1] / math.sqrt(cov_mat_absolute_belowcut[0][0] * cov_mat_absolute_belowcut[1][1]))

    # Plot scatter, and label
    plt.plot(np.absolute(module_misalignments_sample), chi_sq_vals_sample, 'bo')
    plt.ylim(ymax=chi_sq_cut)
    plt.title("Plot of $\chi^2_{red}$ Against Modulus of Module Misalignment, \n for Fits With $\chi^2_{red} < 2$, With " + str(dof) + " D.o.f.")
    plt.ylabel("$\chi^2_{red}$")
    plt.xlabel("Absolute Module Misalignment / mm")

    # Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
    number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_belowcut))))
    shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len([i for i in chi_sq_vals_sample if i < chi_sq_cut ]))))
    correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_absolute_belowcut[0][1] / math.sqrt(cov_mat_absolute_belowcut[0][0] * cov_mat_absolute_belowcut[1][1]))))

    # Get handles and labels from patches, and create legend
    handles, labels = plt.subplot(111).get_legend_handles_labels()
    handles.append(number_patch)
    handles.append(shown_patch)
    handles.append(correlation_patch)
    labels.append(number_patch.get_label())
    labels.append(shown_patch.get_label())
    labels.append(correlation_patch.get_label())
    plt.legend(handles, labels, loc='best', frameon=False)

    # Save figure, and clear
    plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_" +str(dof) + "dof_belowcut_absolute.pdf")
    plt.clf()



#
# All Degrees of Freedom
#
#
# Get sample of values from array. If less points than desired, set sample to array of all points
chi_sq_vals_sample = []
module_misalignments_sample = []
if len(module_misalignments_alldof) >= points_shown:
    module_misalignments_sample = module_misalignments_alldof[:points_shown]
    chi_sq_vals_sample = chi_sq_vals_alldof[:points_shown]
else:
    module_misalignments_sample = module_misalignments_alldof
    chi_sq_vals_sample = chi_sq_vals_alldof
        
# Array of values above and below chi-squared cuts
chi_sq_vals_abovecut = []
chi_sq_vals_belowcut = []
module_misalignments_abovecut = []
module_misalignments_belowcut = []
    
# Populate arrays with values above, below cut
for i in xrange(len(chi_sq_vals_alldof)):
    if chi_sq_vals_alldof[i] > chi_sq_cut:
        chi_sq_vals_abovecut.append(chi_sq_vals_alldof[i])
        module_misalignments_abovecut.append(module_misalignments_alldof[i])
    else:
        chi_sq_vals_belowcut.append(chi_sq_vals_alldof[i])
        module_misalignments_belowcut.append(module_misalignments_alldof[i])
            
# Generate matrices of covariance with and without cut, and for misalignments and modulus of misalignments 
cov_mat = np.cov([module_misalignments_alldof, chi_sq_vals_alldof])
cov_mat_absolute = np.cov([np.absolute(module_misalignments_alldof), chi_sq_vals_alldof])
cov_mat_abovecut = np.cov([module_misalignments_abovecut, chi_sq_vals_abovecut])
cov_mat_absolute_abovecut = np.cov([np.absolute(module_misalignments_abovecut), chi_sq_vals_abovecut])
cov_mat_belowcut = np.cov([module_misalignments_belowcut, chi_sq_vals_belowcut])
cov_mat_absolute_belowcut = np.cov([np.absolute(module_misalignments_belowcut), chi_sq_vals_belowcut])


print ""
print ""
print "All Degrees of Freedom"
print ""
print "No Cut"
print ""
print "Misalignment and Chi-Squared:"
print "Covariance Matrix"
print cov_mat
print ""
print "Correlation: ", str(cov_mat[0][1] / math.sqrt(cov_mat[0][0] * cov_mat[1][1]))
print ""
print ""

# Plot scatter, and label
plt.plot(module_misalignments_sample, chi_sq_vals_sample, 'bo')
plt.title("Plot of $\chi^2_{red}$ Against Module Misalignment, \n for Each Fit, With All D.o.F")
plt.ylabel("$\chi^2_{red}$")
plt.xlabel("Module Misalignment / mm")

# Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_alldof))))
shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len(chi_sq_vals_sample))))
correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat[0][1] / math.sqrt(cov_mat[0][0] * cov_mat[1][1]))))

# Get handles and labels from patches, and create legend
handles, labels = plt.subplot(111).get_legend_handles_labels()
handles.append(number_patch)
handles.append(shown_patch)
handles.append(correlation_patch)
labels.append(number_patch.get_label())
labels.append(shown_patch.get_label())
labels.append(correlation_patch.get_label())
plt.legend(handles, labels, loc='best', frameon=False)

# Save figure, and clear
plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_alldof_nocut.pdf")
plt.clf()


print "Absolute Misalignment and Chi-Squared"
print "Covariance Matrix (Absolute Misalignment)"
print cov_mat_absolute
print ""
print "Correlation: ", str(cov_mat_absolute[0][1] / math.sqrt(cov_mat_absolute[0][0] * cov_mat_absolute[1][1]))

# Plot scatter, and label
plt.plot(np.absolute(module_misalignments_sample), chi_sq_vals_sample, 'bo')
plt.title("Plot of $\chi^2_{red}$ Against Modulus of Module Misalignment, \n for Each Fit, With All D.o.F")
plt.ylabel("$\chi^2_{red}$")
plt.xlabel("Absolute Module Misalignment / mm")

# Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_alldof))))
shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len(chi_sq_vals_sample))))
correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_absolute[0][1] / math.sqrt(cov_mat_absolute[0][0] * cov_mat_absolute[1][1]))))

# Get handles and labels from patches, and create legend
handles, labels = plt.subplot(111).get_legend_handles_labels()
handles.append(number_patch)
handles.append(shown_patch)
handles.append(correlation_patch)
labels.append(number_patch.get_label())
labels.append(shown_patch.get_label())
labels.append(correlation_patch.get_label())
plt.legend(handles, labels, loc='best', frameon=False)

# Save figure, and clear
plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_alldof_nocut_absolute.pdf")
plt.clf()


print ""
print ""
print "Above Cut"
print ""
print "Misalignment and Chi-Squared:"
print "Covariance Matrix"
print cov_mat_abovecut
print ""
print "Correlation: ", str(cov_mat_abovecut[0][1] / math.sqrt(cov_mat_abovecut[0][0] * cov_mat_abovecut[1][1]))
print ""
print ""

# Plot scatter, and label
plt.plot(module_misalignments_sample, chi_sq_vals_sample, 'bo')
plt.ylim(ymin=chi_sq_cut)
plt.title("Plot of $\chi^2_{red}$ Against Module Misalignment, \n for Fits With $\chi^2_{red} > 2$, With All D.o.F")
plt.ylabel("$\chi^2_{red}$")
plt.xlabel("Module Misalignment / mm")

# Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_abovecut))))
shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len([i for i in chi_sq_vals_sample if i > chi_sq_cut ]))))
correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_abovecut[0][1] / math.sqrt(cov_mat_abovecut[0][0] * cov_mat_abovecut[1][1]))))

# Get handles and labels from patches, and create legend
handles, labels = plt.subplot(111).get_legend_handles_labels()
handles.append(number_patch)
handles.append(shown_patch)
handles.append(correlation_patch)
labels.append(number_patch.get_label())
labels.append(shown_patch.get_label())
labels.append(correlation_patch.get_label())
plt.legend(handles, labels, loc='best', frameon=False)

# Save figure, and clear
plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_alldof_abovecut.pdf")
plt.clf()


print "Absolute Misalignment and Chi-Squared"
print "Covariance Matrix (Absolute Misalignment)"
print cov_mat_absolute_abovecut
print ""
print "Correlation: ", str(cov_mat_absolute_abovecut[0][1] / math.sqrt(cov_mat_absolute_abovecut[0][0] * cov_mat_absolute_abovecut[1][1]))

# Plot scatter, and label
plt.plot(np.absolute(module_misalignments_sample), chi_sq_vals_sample, 'bo')
plt.title("Plot of $\chi^2_{red}$ Against Modulus of Module Misalignment, \n for Fits With $\chi^2_{red} > 2$, With All D.o.F")
plt.ylabel("$\chi^2_{red}$")
plt.xlabel("Absolute Module Misalignment / mm")
plt.ylim(ymin=chi_sq_cut)

# Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_abovecut))))
shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len([i for i in chi_sq_vals_sample if i > chi_sq_cut ]))))
correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_absolute_abovecut[0][1] / math.sqrt(cov_mat_absolute_abovecut[0][0] * cov_mat_absolute_abovecut[1][1]))))

# Get handles and labels from patches, and create legend
handles, labels = plt.subplot(111).get_legend_handles_labels()
handles.append(number_patch)
handles.append(shown_patch)
handles.append(correlation_patch)
labels.append(number_patch.get_label())
labels.append(shown_patch.get_label())
labels.append(correlation_patch.get_label())
plt.legend(handles, labels, loc='best', frameon=False)

# Save figure, and clear
plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_alldof_abovecut_absolute.pdf")
plt.clf()


print ""
print ""
print "Below Cut"
print ""
print "Misalignment and Chi-Squared:"
print "Covariance Matrix"
print cov_mat_belowcut
print ""
print "Correlation: ", str(cov_mat_belowcut[0][1] / math.sqrt(cov_mat_belowcut[0][0] * cov_mat_belowcut[1][1]))
print ""
print ""

# Plot scatter, and label
plt.plot(module_misalignments_sample, chi_sq_vals_sample, 'bo')
plt.ylim(ymax=chi_sq_cut)
plt.title("Plot of $\chi^2_{red}$ Against Module Misalignment, \n for Fits With $\chi^2_{red} < 2$, With All D.o.F")
plt.ylabel("$\chi^2_{red}$")
plt.xlabel("Module Misalignment / mm")

# Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_belowcut))))
shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len([i for i in chi_sq_vals_sample if i < chi_sq_cut ]))))
correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_belowcut[0][1] / math.sqrt(cov_mat_belowcut[0][0] * cov_mat_belowcut[1][1]))))

# Get handles and labels from patches, and create legend
handles, labels = plt.subplot(111).get_legend_handles_labels()
handles.append(number_patch)
handles.append(shown_patch)
handles.append(correlation_patch)
labels.append(number_patch.get_label())
labels.append(shown_patch.get_label())
labels.append(correlation_patch.get_label())
plt.legend(handles, labels, loc='best', frameon=False)

# Save figure, and clear
plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_alldof_belowcut.pdf")
plt.clf()


print "Absolute Misalignment and Chi-Squared"
print "Covariance Matrix (Absolute Misalignment)"
print cov_mat_absolute_belowcut
print ""
print "Correlation: ", str(cov_mat_absolute_belowcut[0][1] / math.sqrt(cov_mat_absolute_belowcut[0][0] * cov_mat_absolute_belowcut[1][1]))

# Plot scatter, and label
plt.plot(np.absolute(module_misalignments_sample), chi_sq_vals_sample, 'bo')
plt.ylim(ymax=chi_sq_cut)
plt.title("Plot of $\chi^2_{red}$ Against Modulus of Module Misalignment, \n for Fits With $\chi^2_{red} < 2$, With All D.o.F")
plt.ylabel("$\chi^2_{red}$")
plt.xlabel("Absolute Module Misalignment / mm")

# Create patches showing total data count to find covariance matrix, number of points shown, and correlation found from covariance matrix
number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals_belowcut))))
shown_patch = mpatches.Patch(color='none', label=("Points Shown: " + str(len([i for i in chi_sq_vals_sample if i < chi_sq_cut ]))))
correlation_patch = mpatches.Patch(color='none', label=("Correlation " + str(cov_mat_absolute_belowcut[0][1] / math.sqrt(cov_mat_absolute_belowcut[0][0] * cov_mat_absolute_belowcut[1][1]))))

# Get handles and labels from patches, and create legend
handles, labels = plt.subplot(111).get_legend_handles_labels()
handles.append(number_patch)
handles.append(shown_patch)
handles.append(correlation_patch)
labels.append(number_patch.get_label())
labels.append(shown_patch.get_label())
labels.append(correlation_patch.get_label())
plt.legend(handles, labels, loc='best', frameon=False)

# Save figure, and clear
plt.savefig(output_dir + "chi_corr_scatter_" + straw_status + "_" + str(track_count) + "track_alldof_belowcut_absolute.pdf")
plt.clf()

