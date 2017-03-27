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


track_count = 0
fit_count = 0
fits_attempted = 10000
missing_straws = False
input_dir = ""
output_dir = ""


argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "h:i:t:o:m:f:", ["help", "input_dir=", "track_count=", "output_dir=", "missing_straws=", "fits_attempted="])
except getopt.GetoptError:
    print "tree_reader.py -i <input_dir> -t <track_count> -f <fit_count> -o <output_dir> -m <missing_straws>"
    sys.exit(2)

# Set parameters according to arguments
for opt, arg in opts:

    if opt in ("-h", "--help"):
        print "tree_reader.py -i <input_dir> -t <track_count> -o <output_dir>"
        sys.exit()
    elif opt in ("-i", "--input_dir"):
        input_dir = arg
    elif opt in ("-t", "--track_count"):
        track_count = int(arg)
    elif opt in ("-o", "--output_dir"):
        output_dir = arg
    elif opt in ("-m", "--missing_straws"):
        if arg == "True":
            missing_straws = True
        elif arg == "False":
            missing_straws = False
    elif opt in ("-f", "--fits_attempted"):
        fits_attempted = int(arg)


# String for input file name
straw_status = "missed_straws" if missing_straws else "all_straws"
input_file = input_dir + "nosmear_" + straw_status + "_" + str(track_count) + "_track_" + str(fits_attempted) + "_fit.root"

print input_file

# Get input file
f = ROOT.TFile.Open(input_file, "read")

# Get tree of straw hit events
tree = f.Get("event_tree")

# Lists of module alignments, and dictionary of lists of reduced chi-squared values for different degreees of freedom
module_true_alignments = []
module_fitted_alignments = []

# Get number of fits in tree, by examining last entry in tree 
entry_count = tree.GetEntries()
tree.GetEntry(entry_count - 1)
fit_count = tree.fitNum + 1



print "Total Fits:", str(fit_count)
#fit_cut_tree = 0        

for i in xrange(fit_count):
    
    if (i % 100) == 0:
        print i
        #memtrack.print_diff()

    # Get entries with correct fit number in tree
    fit_cut_tree = tree.CopyTree('fitNum==' + str(i))

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


    fit_cut_tree.GetEntry(0)

    # Append module alignments to lists
    module_true_alignments.append(fit_cut_tree.trueAlignment)
    module_fitted_alignments.append(fit_cut_tree.fittedAlignment)

    fit_cut_tree.IsA().Destructor(fit_cut_tree) # To avoid stupid bloody memory leak

# Calculate array of module misalignments
module_misalignments = np.array(module_fitted_alignments) - np.array(module_true_alignments)

# Calculate mean, standard deviation, and number of data points in array of module misalignments
align_mean = np.mean(module_misalignments)
align_std_dev = np.std(module_misalignments)
align_data_count = len(module_misalignments)


# Plot histogram of alignment errors
plt.hist(module_misalignments, bins=20, color='blue')
plt.title("Plot of Module Misalignment Distances, With All D.o.F.")
plt.xlabel("Misalignment / mm")
plt.ylabel("Frequency")

# Entries to put in legend, with relevent statistical quantities
number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(align_data_count)))
mean_patch = mpatches.Patch(color='none', label=("$\mu$: " + str(align_mean)))
err_mean_patch = mpatches.Patch(color='none', label=("$SE(\mu)$: " + str(align_std_dev / math.sqrt(align_data_count))))
std_patch = mpatches.Patch(color='none', label=("$\sigma$: " + str(align_std_dev)))
err_std_patch = mpatches.Patch(color='none', label=("$SE(\sigma)$: " + str(align_std_dev / math.sqrt((2 * align_data_count) - 2))))

# Get legend handles and labels, and plot legend
align_handles = [number_patch, mean_patch, err_mean_patch, std_patch, err_std_patch]
align_labels = [number_patch.get_label(), mean_patch.get_label(), err_mean_patch.get_label(), std_patch.get_label(), err_std_patch.get_label()]
plt.legend(align_handles, align_labels, loc='best', frameon=False)

# Save plot
plt.savefig(output_dir + "nosmear_" + straw_status + "_" + str(track_count) + "t_align.pdf")
plt.clf()


f.Close()

