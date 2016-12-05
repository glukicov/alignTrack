import toytracebackthreemod as ttb
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
chi_sq_cut = 2
missing_straws = False
input_dir = ""
output_dir = ""


argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "h:i:t:c:o:m:f:", ["help", "input_dir=", "track_count=", "chi_sq_cut=", "output_dir=", "missing_straws=", "fits_attempted="])
except getopt.GetoptError:
    print "tree_reader.py -i <input_dir> -t <track_count> -f <fit_count> -c <chi_sq_cut> -o <output_dir> -m <missing_straws>"
    sys.exit(2)

# Set parameters according to arguments
for opt, arg in opts:

    if opt in ("-h", "--help"):
        print "tree_reader.py -i <input_dir> -t <track_count> -c <chi_sq_cut> -o <output_dir>"
        sys.exit()
    elif opt in ("-i", "--input_dir"):
        input_dir = arg
    elif opt in ("-t", "--track_count"):
        track_count = int(arg)
    elif opt in ("-c", "--chi_sq_cut"):
        chi_sq_cut = float(arg)
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
input_file = input_dir + "smear_" + straw_status + "_" + str(track_count) + "_track_" + str(fits_attempted) + "_fit.root"

print input_file

# Get input file
f = ROOT.TFile.Open(input_file, "read")

# Get tree of straw hit events
tree = f.Get("event_tree")

# Dictionary of lists of reduced chi-squared values for different degrees of freedom
chi_sq_vals = dict()

# 2D Histogram of misalignments for each module
hMisalignments = ROOT.TH2D("hMisalignments", "hMisalignments", 100, -2.5, 2.5, 100, -2.5, 2.5)

# Get number of fits in tree, by examining last entry in tree 
entry_count = tree.GetEntries()
tree.GetEntry(entry_count - 1)
fit_count = tree.fitNum + 1


print "Total Fits:", str(fit_count)
        
for i in xrange(fit_count):
    
    if (i % 100) == 0:
        print i

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
    
    # Create array of hit distance uncertainties, taken to be equal to hit resolution of detector
    sigmas = np.array([ttb.hit_resolution for n in xrange(fit_cut_tree_entry_count)])

    # Calculate reduced value of chi-squared from hit distances
    chi_sq_red = np.sum(((np.array(true_hit_distances) - np.array(fitted_hit_distances)) / sigmas)**2) / (fit_cut_tree_entry_count - ((2 * track_count) + 2)) 

    # Check if chi-squared below cut
    if (chi_sq_red < chi_sq_cut):
        fit_cut_tree.GetEntry(0)

        # Add entry in 2D histogram for misalignments
        hMisalignments.Fill((fit_cut_tree.fittedMod1Alignment - fit_cut_tree.trueMod1Alignment), (fit_cut_tree.fittedMod2Alignment - fit_cut_tree.trueMod2Alignment))

        # Add chi-squared value to dictionary
        if fit_cut_tree_entry_count - ((2 * track_count) + 2) in chi_sq_vals:
            chi_sq_vals[fit_cut_tree_entry_count - ((2 * track_count) + 2)].append(chi_sq_red)
        else:
            chi_sq_vals[fit_cut_tree_entry_count - ((2 * track_count) + 2)] = [chi_sq_red]


# Draw 2D histogram, sorting out cosmetics, then print
cMisalignments = ROOT.TCanvas("cMisalignments", "cMisalignments", 2000, 1500)
ROOT.gStyle.SetPalette(51)
ROOT.gStyle.SetOptStat(2210)
hMisalignments.GetXaxis().SetTitle("Module 1 Misalignment / mm")
hMisalignments.GetYaxis().SetTitle("Module 2 Misalignment / mm")
hMisalignments.SetTitle("Misalignments in Each Module, \n for Fits With " + str(track_count) + " Tracks (chi^2 < " + str(chi_sq_cut) + ")")
hMisalignments.Draw("colz")

ROOT.gStyle.SetStatX(0.9)
ROOT.gStyle.SetStatY(0.9)
ROOT.gStyle.SetStatW(0.15)
ROOT.gStyle.SetStatH(0.2)

cMisalignments.SetCanvasSize(2000, 1500)
cMisalignments.Print(output_dir + straw_status + "_" + str(track_count) + "t_align.pdf")

# Loop over degrees of freedom of recorded chi-squared values
for dof in chi_sq_vals:
    
    # Draw histogram of chi-squared values, getting bin contents and positions for fitting to ideal distributions
    chi_bin_conts, chi_bin_edges, chi_patches = plt.hist(chi_sq_vals[dof], bins=50, color='blue', label="Calculated $\chi^2_{red}$")
    chi_bin_width = (chi_bin_edges[-1] - chi_bin_edges[0]) / len(chi_bin_conts)

    # Plot ideal chi-squared distribution
    x_chi_plot = np.linspace(0, 2, 100)
    plt.plot(x_chi_plot, len(chi_sq_vals[dof]) * chi_bin_width * stats.chi2.pdf(x_chi_plot, dof, scale=1.0/dof), "-r", label=("Chi-Squared Dist (" + str(dof) + " dof)"))

    # Create arrays of number of fits with chi-squared values within given bin, and number expected in given bin from chi-squared dist.
    chi_func_bins_unfiltered = np.array([len(chi_sq_vals[dof]) * (stats.chi2.cdf(chi_bin_edges[i+1], dof) - stats.chi2.cdf(chi_bin_edges[i], dof)) for i in xrange(len(chi_bin_conts))])
    chi_obs_bins_unfiltered = np.array(chi_bin_conts)
    

    # Filter out entries in arrays corresponding to bins with zero entries
    bin_count = len(chi_func_bins_unfiltered)
    chi_func_bins = np.array([chi_func_bins_unfiltered[x] for x in xrange(bin_count) if chi_obs_bins_unfiltered[x] > 0])
    chi_obs_bins = np.array([chi_obs_bins_unfiltered[x] for x in xrange(bin_count) if chi_obs_bins_unfiltered[x] > 0])

    chi_sigmas = np.sqrt(chi_obs_bins) # Calculate uncertainty of observed entries in histogram bins. (Assume Poisson stats for each bin.)

    # Calculate chi-squared between observed fit chi-squared values, and ideal distribution
    chi_squared_dist_chi_squared = np.sum(((chi_func_bins - chi_obs_bins) / chi_sigmas)**2)

    # Patches to show chi-squared, number of entries in legend
    number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_sq_vals[dof]))))
    chi_squared_patch = mpatches.Patch(color='none', label=("$\chi^2_{red}$: " + str(chi_squared_dist_chi_squared / len(chi_obs_bins))))
        
    # Sort out handles and labels, then create legend
    chi_handles, chi_labels = plt.subplot(111).get_legend_handles_labels()
    chi_handles.append(number_patch)
    chi_handles.append(chi_squared_patch)
    chi_labels.append(number_patch.get_label())
    chi_labels.append(chi_squared_patch.get_label())
    plt.legend(chi_handles, chi_labels, loc='best', frameon=False)

    # Sort out titles
    plt.xlabel("Fit $\chi^2_{red}$")
    plt.ylabel("Frequency")
    plt.title("Plot of $\chi^2$ for Fits With $\chi^2_{red} < " + str(chi_sq_cut) + "$, With " + str(dof) + " D.o.F., Compared to Ideal Distribution")

    # Save plot
    plt.savefig(output_dir + straw_status + "_" + str(track_count) + "t_chi" + str(dof) + ".pdf")
    plt.clf()

f.Close()

