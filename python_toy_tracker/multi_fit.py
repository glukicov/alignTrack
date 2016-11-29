import toytraceback as ttb
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import random
import sys
import math
import scipy.stats as stats
import time
import os.path
import getopt
from scipy.optimize import curve_fit 

# Script to plot alignment residuals and chi-squared distributions for multiple requested track counts, and a requested number
# of fits. Intended for use of batch system.

# Get start CPU time, and wall time (in seconds), to ensure script doesn't overrun batch system limits
start_wall_time = time.time()
start_cpu_time = time.clock()

track_counts = [5, 10, 15, 20, 30, 40, 50] # Number of tracks to fit
module_alignment = 1.5 # Alignment of first module (in mm)
fit_count = 10 # Number of times to run fit

# Whether to include smearing of residuals, and finite straw size
detector_smear = True
detector_missing_hits = False

# Time limits for script on batch system (in s), with time buffer to ensure partial results are processed. 
wall_time_limit = 345600
cpu_time_limit = 259200
time_buffer = 100

# Dictionary of CPU and wall times for different UCL Batch queues
queue_time_dict = {"short": [3600, 3600], "medium": [28800, 86400], "long": [259200, 345600]}

# Get arguments and options from console
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "hf:s:m:q:a:o:", ["fit_count=", "smearing=", "missing_hits=", "batch_queue=", "module_alignment=", "output_dir="])
except  getopt.GetoptError:
    print "multi_fit.py -f <fit_count> -s <smearing> -m <missing_hits> -q <batch_queue> -a <module_alignment> -o <output_dir>"
    sys.exit(2)

# Set parameters according to arguments
for opt, arg in opts:
    if opt == "-h":
        print "multi_fit.py -f <fit_count> -s <smearing> -m <missing_hits> -q <batch_queue> -a <module_alignment> -o <output_dir>"
        sys.exit()
    elif opt in ("-f", "--fit_count"):
        fit_count = int(arg)
    elif opt in ("-s", "--smearing"):
        if (arg == "True"):
            detector_smear = True
        elif (arg == "False"):
            detector_smear = False
        else:
            print "-s, --smearing: Must use 'True' or 'False'"
            sys.exit()
    elif opt in ("-m", "--missing_hits"):
        if (arg == "True"):
            detector_missing_hits = True
        elif (arg == "False"):
            detector_missing_hits = False
        else:
            print "-m, --missing_hits: Must use 'True' or 'False'"
            sys.exit()
    elif opt in ("-q", "--batch_queue"):
        if (arg in ("short", "medium", "long")):
            cpu_time_limit = queue_time_dict[arg][0]
            wall_time_limit = queue_time_dict[arg][1]
        else:
            print "-q, --batch_queue: Must use 'short', 'medium', or 'long'"
            sys.exit()
    elif opt in ("-a", "-module_alignment"):
        module_alignment = float(arg)
    

# Loop across track counts to be examined by script
for track_count in track_counts:

    # Contain calculated chi^2 values, and alignment residuals. Use dictionary with lists of values for each number of degrees of freedom in fit. 
    alignment_errors = dict()
    residual_means = dict()
    chi_squared_fit_vals = dict()
    chi_squared_res_vals = dict()

    # Maximum reduced chi-squared value, to exclude poor fits. Fits with value above this will be discarded. 
    reduced_chi_squared_cut = 2
            
    print "Carrying out", fit_count, "fits, with", track_count, "tracks:"  
        
    # Carry out requested number of fits
    for k in xrange(fit_count):
    
        print k
        
        # FITTING

        # Create instance of detector, and set alignment to desired value
        true_detector = ttb.Detector(smearing=detector_smear, finite_straws=detector_missing_hits)
        true_detector.set_module_x_align(1, module_alignment)

        # # Get wire coordinates
        # x_true_wires = true_detector.get_wires_x()
        # y_true_wires = true_detector.get_wires_y()

        # Coordinates at beginning, end of track
        x_bottom = [random.uniform(ttb.track_boundary_x_lower, ttb.track_boundary_x_upper) for i in xrange(track_count)]
        x_top = [random.uniform(ttb.track_boundary_x_lower, ttb.track_boundary_x_upper) for i in xrange(track_count)]
        y_bottom = [ttb.track_boundary_y_lower for i in xrange(track_count)]
        y_top = [ttb.track_boundary_y_upper for i in xrange(track_count)]
    
        
        # Calculate track gradient and x-intercept from track coordinates
        gradient = [((x_top[i] - x_bottom[i]) / (y_top[i] - y_bottom[i])) for i in xrange(track_count)]
        intercept = [x_bottom[i] - (gradient[i] * y_bottom[i]) for i in xrange(track_count)]
        
        # Set up track, and get pos of beginning, end points
        true_track = [ttb.Track(gradient[i], intercept[i]) for i in xrange(track_count)]
        # x_true_track = [true_track[i].get_x_points() for i in xrange(track_count)]
        # y_true_track = [true_track[i].get_y_points() for i in xrange(track_count)]
        
        
        # Get wire hits in each layer, then assign these to detector
        events = []
        for j in xrange(track_count): 
            wire_hits = ttb.closest_hit_wires(true_detector, true_track[j])
            event = ttb.DetectorHitEvent(true_track[j], wire_hits, j)
            events.append(event)
            

        # New detector object for fitting
        fitting_detector = ttb.Detector(smearing=detector_smear, finite_straws=detector_missing_hits)
        fitting_detector.set_events(events)

        # Get x-displacement of all closest approached wires
        hit_rads = []
        for event in events:
            for wire_hit in event.wire_hits:
                hit_rads.append(wire_hit.get_hit_dist())
                
        # Get list of keys for all struck wires, by all tracks
        wire_keys = fitting_detector.get_hit_keys()

        guess = [0.0 for i in xrange((2 * track_count) + 1)] # Initial guess for fitted track gradient and intercept, and module alignment (spurious convergence without this)
    
    
        fit_sigmas = [ttb.hit_resolution]*len(hit_rads) # Uncertainties on hit radii, from detector hit resolution

        # Finds values for track gradient and intercept, by fitting to wire hit x-displacements
        popt, pcov = curve_fit(fitting_detector.get_hit_radius, wire_keys, hit_rads, p0=guess, sigma=fit_sigmas)

        # Calculate residuals and chi-squared for residuals, to test goodness of fit.
        residuals = np.array(hit_rads) - np.array(fitting_detector.get_hit_radius(wire_keys, popt)) / 1.0
        deg_freedom = len(residuals) - len(popt)
        chi_squared = np.sum((residuals / fit_sigmas)**2)


        # FITTING GAUSSIAN TO HIT RESIDUALS

        # Filter large residuals for fitting of gaussian function
        filtered_residuals = [x for x in residuals if abs(x) < (10 * ttb.hit_resolution)]

        # Get mean and std-dev of fitted gaussian, then array of x-values for plotting
        (mu_res, sigma_res) = stats.norm.fit(filtered_residuals)
        gaus_x_res = np.linspace(mu_res - 4 * sigma_res, mu_res + 4 * sigma_res, 500)

        # Plot histogram of residuals, with large ones filtered out
        res_obs_bin_contents_unfiltered, res_obs_bin_edges, _ = plt.hist(filtered_residuals, bins=20)

        # Get expected number of residuals in each bin, from fitted gaussian
        bin_count = len(res_obs_bin_contents_unfiltered)
        res_func_bin_contents_unfiltered = [len(filtered_residuals) * (stats.norm.cdf(res_obs_bin_edges[i+1], scale=(sigma_res), loc=mu_res) - stats.norm.cdf(res_obs_bin_edges[i], scale=(sigma_res), loc=mu_res)) for i in xrange(bin_count)]

        # Filter out bins containing zero observed residuals, to exclude from chi-squared calculation
        res_fun_bin_contents = np.array([res_func_bin_contents_unfiltered[i] for i in xrange(bin_count) if res_obs_bin_contents_unfiltered[i] > 0])
        res_obs_bin_contents = np.array([res_obs_bin_contents_unfiltered[i] for i in xrange(bin_count) if res_obs_bin_contents_unfiltered[i] > 0])

        res_sigmas = np.sqrt(res_fun_bin_contents) # Uncertainty in observed bins (Poisson)

        # Calculate chi-squared for histogram of residuals and gaussian fit
        chi_squared_fitted_gauss = np.sum(((res_fun_bin_contents - res_obs_bin_contents) / res_sigmas)**2) / bin_count
        
        # Clear plot of histogram
        plt.clf()


        # RECORDING RESULTS

        # Get alignment residual, and chi-squared value if reduced chi-squared less than 10
        # (Should filter out spurious fits)
        if (chi_squared < reduced_chi_squared_cut * deg_freedom):

            # Check if list of values exists for this number of degrees of freedom. Create or append as necessary
            if (deg_freedom in alignment_errors):
                chi_squared_fit_vals[deg_freedom].append(chi_squared) 
                chi_squared_res_vals[deg_freedom].append(chi_squared_fitted_gauss)
                alignment_errors[deg_freedom].append(popt[-1] - module_alignment)
                residual_means[deg_freedom].append(mu_res)
            else:
                chi_squared_fit_vals[deg_freedom] = [chi_squared]
                alignment_errors[deg_freedom] = [popt[-1] - module_alignment]
                chi_squared_res_vals[deg_freedom] = [chi_squared_fitted_gauss] 
                residual_means[deg_freedom] = [mu_res]

        # If close to time limit, break out of loop of fits, so partial results can be processed
        if (time.time() - start_wall_time) > (wall_time_limit - time_buffer):
            fit_count = k + 1
            break
        if (time.clock() - start_cpu_time) > (cpu_time_limit - time_buffer):
            fit_count = k + 1
            break


    # Set name of directory to save histograms in, with information on number of fits and tracks, whether finite straw sizes are used, whether smearing is implemented  
    smear_string = "smear" if detector_smear else "no_smear"
    straw_string = "missed_straws" if detector_missing_hits else "all_straws"
    directory = smear_string + "_" + straw_string + "_" + str(track_count) + "_track_" + str(fit_count) + "_fit/"

    # Create directory if doesn't already exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    # List to contain alignment errors and residual means for all degrees of freedom.
    align_errors_all_dof = []
    residual_means_all_dof = []
    chi_squared_fit_all_dof = []
    chi_squared_res_all_dof = []

    print ""
    print ""
    # Loop over numbers of degrees of freedom recorded
    for deg_freedom in alignment_errors:

        # ALIGNMENT ERRORS

        # append alignment errors with this number of degrees of freedom to list containing errors with all dof.
        align_errors_all_dof.extend(alignment_errors[deg_freedom])

        # Calculate mean, standard deviation, and number of data points for alignment errors
        align_mean = np.mean(np.array(alignment_errors[deg_freedom]))
        align_std_dev = np.std(np.array(alignment_errors[deg_freedom]))
        align_data_count = len(alignment_errors[deg_freedom])

        print "*******************************************"
        print "Fits With", deg_freedom, "Degrees of Freedom"
        print "*******************************************"
        print ""
        print "Mean Alignment Error:", align_mean, "mm"
        print "(Standard Error):", align_std_dev / math.sqrt(align_data_count), "mm" 
        print "Standard Deviation:", align_std_dev, "mm"
        print "(Standard Error):", align_std_dev / math.sqrt((2 * align_data_count) - 2), "mm"

        # Plot histogram of alignment errors
        plt.hist(alignment_errors[deg_freedom], bins=20, color='blue')
        plt.title("Plot of Alignment Fit Errors for Fits With $\chi^2_{red} < " + str(reduced_chi_squared_cut) + "$, With " + str(deg_freedom) + " D.o.F.")
        plt.xlabel("Alignment Error Value / mm")
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

        # Set appropriate file name, and save plot
        file_name_prefixes = directory + str(deg_freedom) + "_dof_"
        plt.savefig(file_name_prefixes + "align.pdf", format = "pdf")
        plt.clf()
    

        # RESIDUAL MEANS

        residual_means_all_dof.extend(residual_means[deg_freedom])

        # Calculate mean, standard deviation, and number of data points for hit residual means
        res_mean_mean = np.mean(np.array(residual_means[deg_freedom]))
        res_mean_std_dev = np.std(np.array(residual_means[deg_freedom]))
        res_mean_data_count = len(residual_means[deg_freedom])

        print ""
        print "Mean of Means of Gaussians Fitted to Hit Residuals:", res_mean_mean, "mm"
        print "(Standard Error):", res_mean_std_dev / math.sqrt(res_mean_data_count), "mm" 
        print "Standard Deviation:", res_mean_std_dev, "mm"
        print "(Standard Error):", res_mean_std_dev / math.sqrt((2 * res_mean_data_count) - 2), "mm"

        # Plot histogram of hit residual means
        plt.hist(residual_means[deg_freedom], bins=20, color='blue')
        plt.title("Plot of Means of Gaussians Fitted to Hit Residuals for \n Event Fits With $\chi^2_{red} < " + str(reduced_chi_squared_cut) + "$, With " + str(deg_freedom) + " D.o.F.")
        plt.xlabel("Hit Residual Fitted Gaussian Mean / mm")
        plt.ylabel("Frequency")

        # Entries to put in legend, with relevent statistical quantities
        number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(res_mean_data_count)))
        mean_patch = mpatches.Patch(color='none', label=("$\mu$: " + str(res_mean_mean)))
        err_mean_patch = mpatches.Patch(color='none', label=("$SE(\mu)$: " + str(res_mean_std_dev / math.sqrt(res_mean_data_count))))
        std_patch = mpatches.Patch(color='none', label=("$\sigma$: " + str(res_mean_std_dev)))
        err_std_patch = mpatches.Patch(color='none', label=("$SE(\sigma)$: " + str(res_mean_std_dev / math.sqrt((2 * res_mean_data_count) - 2))))

        # Get legend handles and labels, and plot legend
        res_mean_handles = [number_patch, mean_patch, err_mean_patch, std_patch, err_std_patch]
        res_mean_labels = [number_patch.get_label(), mean_patch.get_label(), err_mean_patch.get_label(), std_patch.get_label(), err_std_patch.get_label()]
        plt.legend(res_mean_handles, res_mean_labels, loc='best', frameon=False)

        # Set appropriate file name, and save plot
        file_name_prefixes = directory + str(deg_freedom) + "_dof_"
        plt.savefig(file_name_prefixes + "res_mean.pdf", format = "pdf")
        plt.clf()


        # CHI-SQUARED VALUES

        chi_squared_fit_all_dof.extend([x / deg_freedom for x in chi_squared_fit_vals[deg_freedom]])

        # Plot histogram of chi-squared, then get bin contents, and edges
        x_chi_plot = np.linspace(stats.chi2.ppf(0.01, deg_freedom), stats.chi2.ppf(0.99, deg_freedom), 100)
        chi_bin_conts, chi_bin_edges, chi_patches = plt.hist(chi_squared_fit_vals[deg_freedom], bins=20, label='Calculated Chi-Squared', color='blue')

        # Get width of each histogram bin
        chi_bin_width = (chi_bin_edges[-1] - chi_bin_edges[0]) / len(chi_bin_conts)

        # Plot chi-squared distribution for given number of degrees of freedom in fit.
        chi_dist_plot = plt.plot(x_chi_plot, len(chi_squared_fit_vals[deg_freedom]) * chi_bin_width * stats.chi2.pdf(x_chi_plot, deg_freedom), 'r-', label=("Chi-Squared Dist (" + str(deg_freedom) + " dof)"))
        plt.xlabel("Chi-Squared")
        plt.ylabel("Frequency")
        plt.title("Plot of $\chi^2$ for Fits With $\chi^2_{red} < " + str(reduced_chi_squared_cut) + "$, With " + str(deg_freedom) + " D.o.F., Compared to Ideal Distribution")

    
        # Create arrays of number of fits with chi-squared values within given bin, and number expected in given bin from chi-squared dist.
        chi_func_bins_unfiltered = np.array([len(chi_squared_fit_vals[deg_freedom]) * (stats.chi2.cdf(chi_bin_edges[i+1], deg_freedom) - stats.chi2.cdf(chi_bin_edges[i], deg_freedom)) for i in xrange(len(chi_bin_conts))])
        chi_obs_bins_unfiltered = np.array(chi_bin_conts)


        # Filter out entries in arrays corresponding to bins with zero entries
        bin_count = len(chi_func_bins_unfiltered)
        chi_func_bins = np.array([chi_func_bins_unfiltered[x] for x in xrange(bin_count) if chi_obs_bins_unfiltered[x] > 0])
        chi_obs_bins = np.array([chi_obs_bins_unfiltered[x] for x in xrange(bin_count) if chi_obs_bins_unfiltered[x] > 0])

        chi_sigmas = np.sqrt(chi_obs_bins) # Calculate uncertainty of observed entries in histogram bins. (Assume Poisson stats for each bin.)

        # Calculate chi-squared between observed fit chi-squared values, and ideal distribution
        chi_squared_dist_chi_squared = np.sum(((chi_func_bins - chi_obs_bins) / chi_sigmas)**2)

        print ""
        print "Number of bins used for test:", len(chi_obs_bins)
        print "Chi-squared of measured chi-squared values to chi-squared distribution:", chi_squared_dist_chi_squared
        print "P(chi^2):", 1 - stats.chi2.cdf(chi_squared_dist_chi_squared, len(chi_obs_bins))
        print ""

        # Patches to show chi-squared, number of entries in legend
        number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(len(chi_squared_fit_vals[deg_freedom]))))
        chi_squared_patch = mpatches.Patch(color='none', label=("$\chi^2_{red}$: " + str(chi_squared_dist_chi_squared / len(chi_obs_bins))))
        
        # Sort out handles and labels, then create legend
        chi_handles, chi_labels = plt.subplot(111).get_legend_handles_labels()
        chi_handles.append(number_patch)
        chi_handles.append(chi_squared_patch)
        chi_labels.append(number_patch.get_label())
        chi_labels.append(chi_squared_patch.get_label())
        plt.legend(chi_handles, chi_labels, loc='best', frameon=False)

        # Save plot
        plt.savefig(file_name_prefixes + "chi.pdf", format = "pdf")
        plt.clf()

        
        # SCATTER PLOT - CHI-SQUARED EVENT, RESIDUAL GAUSSIAN FIT

        # Append chi-squared values of gaussian fits to residuals for this d.o.f to list for all d.o.f.
        chi_squared_res_all_dof.extend(chi_squared_res_vals[deg_freedom])
        
        # Create scatter plot
        plt.scatter([x / deg_freedom for x in chi_squared_fit_vals[deg_freedom]], chi_squared_res_vals[deg_freedom])
        plt.xlim(0, 2.5)
        plt.ylim(0, 2.5)

        # Label scatter plot
        plt.xlabel("Event Fit $\chi^2_{red}$")
        plt.ylabel("Residuals Gaussian Fit $\chi^2_{red}$")
        plt.title("Scatter Plot of $\chi^2_{red}$ Values for Gaussian Fits to Residuals, \n Against $\chi^2_{red}$ for Event, Using Event Fits With $\chi^2_{red}$ < " + str(reduced_chi_squared_cut))

        # Save plot, then clear
        plt.savefig(file_name_prefixes + "scatter.pdf", color='blue', format = "pdf", )
        plt.clf()


    # Break out of loop of track numbers to fit, to process partial results, if near end of time limit
    if (time.time() - start_wall_time) > (wall_time_limit - time_buffer):
        break
    if (time.clock() - start_cpu_time) > (cpu_time_limit - time_buffer):
        break


    # ALIGNMENT ERRORS - ALL D.O.F.

    # Calculate mean, standard deviation, number of recorded alignment errors
    align_mean = np.mean(np.array(align_errors_all_dof))
    align_std_dev = np.std(np.array(align_errors_all_dof))
    align_data_count = len(align_errors_all_dof)

    print "*******************************************"
    print "Fits With All Degrees of Freedom"
    print "*******************************************"
    print ""
    print "Mean Alignment Error:", align_mean, "mm"
    print "(Standard Error):", align_std_dev / math.sqrt(align_data_count), "mm" 
    print "Standard Deviation:", align_std_dev, "mm"
    print "(Standard Error):", align_std_dev / math.sqrt((2 * align_data_count) - 2), "mm"
    print ""
    
    # Plot histogram of residuals
    plt.hist(align_errors_all_dof, bins=20, color='blue')
    plt.title("Plot of Alignment Fit Errors for Fits With $\chi^2_{red} < " + str(reduced_chi_squared_cut) + "$, With All D.o.F.")
    plt.xlabel("Alignment Error Value / mm")
    plt.ylabel("Frequency")
    
    # To show relevent statistical quantities in legend
    number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(align_data_count)))
    mean_patch = mpatches.Patch(color='none', label=("$\mu$: " + str(align_mean)))
    err_mean_patch = mpatches.Patch(color='none', label=("$SE(\mu)$: " + str(align_std_dev / math.sqrt(align_data_count))))
    std_patch = mpatches.Patch(color='none', label=("$\sigma$: " + str(align_std_dev)))
    err_std_patch = mpatches.Patch(color='none', label=("$SE(\sigma)$: " + str(align_std_dev / math.sqrt((2 * align_data_count) - 2))))
    
    # Get legend handles and labels, and show legend
    align_handles = [number_patch, mean_patch, err_mean_patch, std_patch, err_std_patch]
    align_labels = [number_patch.get_label(), mean_patch.get_label(), err_mean_patch.get_label(), std_patch.get_label(), err_std_patch.get_label()]
    plt.legend(align_handles, align_labels, loc='best', frameon=False)
    
    # Save plot
    file_name_prefixes = directory + "all_dof_"
    plt.savefig(file_name_prefixes + "align.pdf", format = "pdf")
    plt.clf()


    # RESIDUAL MEANS

    # Calculate mean, standard deviation, and number of data points for hit residual means
    res_mean_mean = np.mean(np.array(residual_means_all_dof))
    res_mean_std_dev = np.std(np.array(residual_means_all_dof))
    res_mean_data_count = len(residual_means_all_dof)

    print "Mean of Means of Gaussians Fitted to Hit Residuals:", res_mean_mean, "mm"
    print "(Standard Error):", res_mean_std_dev / math.sqrt(res_mean_data_count), "mm" 
    print "Standard Deviation:", res_mean_std_dev, "mm"
    print "(Standard Error):", res_mean_std_dev / math.sqrt((2 * res_mean_data_count) - 2), "mm"
    print ""
    print ""

    # Plot histogram of hit residual means
    plt.hist(residual_means_all_dof, bins=20, color='blue')
    plt.title("Plot of Means of Gaussians Fitted to Hit Residuals for \n Event Fits With $\chi^2_{red} < " + str(reduced_chi_squared_cut) + "$, With All D.o.F.")
    plt.xlabel("Hit Residual Fitted Gaussian Mean / mm")
    plt.ylabel("Frequency")

    # Entries to put in legend, with relevent statistical quantities
    number_patch = mpatches.Patch(color='none', label=("Data Count: " + str(res_mean_data_count)))
    mean_patch = mpatches.Patch(color='none', label=("$\mu$: " + str(res_mean_mean)))
    err_mean_patch = mpatches.Patch(color='none', label=("$SE(\mu)$: " + str(res_mean_std_dev / math.sqrt(res_mean_data_count))))
    std_patch = mpatches.Patch(color='none', label=("$\sigma$: " + str(res_mean_std_dev)))
    err_std_patch = mpatches.Patch(color='none', label=("$SE(\sigma)$: " + str(res_mean_std_dev / math.sqrt((2 * res_mean_data_count) - 2))))

    # Get legend handles and labels, and plot legend
    res_mean_handles = [number_patch, mean_patch, err_mean_patch, std_patch, err_std_patch]
    res_mean_labels = [number_patch.get_label(), mean_patch.get_label(), err_mean_patch.get_label(), std_patch.get_label(), err_std_patch.get_label()]
    plt.legend(res_mean_handles, res_mean_labels, loc='best', frameon=False)

    # Set appropriate file name, and save plot
    file_name_prefixes = directory + "all_dof_"
    plt.savefig(file_name_prefixes + "res_mean.pdf", format = "pdf")
    plt.clf()

    # SCATTER PLOT - CHI-SQUARED EVENT, RESIDUAL GAUSSIAN FIT
        
    # Create scatter plot
    plt.scatter(chi_squared_fit_all_dof, chi_squared_res_all_dof)
    plt.xlim(0, 2.5)
    plt.ylim(0, 2.5)

    # Label scatter plot
    plt.xlabel("Event Fit $\chi^2_{red}$")
    plt.ylabel("Residuals Gaussian Fit $\chi^2_{red}$")
    plt.title("Scatter Plot of $\chi^2_{red}$ Values for Gaussian Fits to Residuals, \n Against $\chi^2_{red}$ for Event, Using Event Fits With $\chi^2_{red}$ < " + str(reduced_chi_squared_cut))

    # Save plot, then clear
    plt.savefig(file_name_prefixes + "scatter.pdf", color='blue', format = "pdf", )
    plt.clf()




# Get end time of script, then output    
end_wall_time = time.time()
end_cpu_time = time.clock()
print ""
print "Completed in " + str(end_wall_time - start_wall_time) + "s."
print "CPU Time Elapsed: " + str(end_cpu_time - start_cpu_time) + "s."
print ""
