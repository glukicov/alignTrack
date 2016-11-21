import toytraceback as ttb
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
import math
import scipy.stats as stats
from scipy.optimize import curve_fit 

track_count = 2 # Number of tracks to fit
module_alignment = 1.5 # Set up detector, and change alignment of first module
fit_count = 500 # Number of times to run fit

deg_freedom = 0 # Number of degrees of freedom in fit

alignment_residuals = []
chi_squared_vals = []

print "Carrying out", fit_count, "fits, with", track_count, "tracks:"  

for k in xrange(fit_count):
    true_detector = ttb.Detector(smearing=True)
    true_detector.set_module_x_align(1, module_alignment)

    # Get wire coordinates
    x_true_wires = true_detector.get_wires_x()
    y_true_wires = true_detector.get_wires_y()

    # Coordinates at beginning, end of track
    x_bottom = [random.uniform(ttb.track_boundary_x_lower, ttb.track_boundary_x_upper) for i in xrange(track_count)]
    x_top = [random.uniform(ttb.track_boundary_x_lower, ttb.track_boundary_x_upper) for i in xrange(track_count)]
    y_bottom = [ttb.track_boundary_y_lower for i in xrange(track_count)]
    y_top = [ttb.track_boundary_y_upper for i in xrange(track_count)]

    # # Get boundaries on intercepts and gradients for fit, with large boundary on module alignment distance
    # lower_fit_bounds = [0]*((2 * track_count) + 1)
    # upper_fit_bounds = [0]*((2 * track_count) + 1)
    # for i in xrange(len(lower_fit_bounds)):
    #     if (i % 2 == 0):
    #         lower_fit_bounds[i] = (ttb.track_boundary_x_lower - ttb.track_boundary_x_upper) / (ttb.track_boundary_y_upper - ttb.track_boundary_y_lower)
    #         upper_fit_bounds[i] = (ttb.track_boundary_x_upper - ttb.track_boundary_x_lower) / (ttb.track_boundary_y_upper - ttb.track_boundary_y_lower)
    #     else:
    #         lower_fit_bounds[i] = ttb.track_boundary_x_lower
    #         upper_fit_bounds[i] = ttb.track_boundary_x_upper
    
    # lower_fit_bounds[-1] = -np.inf
    # upper_fit_bounds[-1] = np.inf
    
    # fit_bounds = (lower_fit_bounds, upper_fit_bounds)
    # print fit_bounds
    
    
    
    # Calculate track gradient and x-intercept from track coordinates
    gradient = [((x_top[i] - x_bottom[i]) / (y_top[i] - y_bottom[i])) for i in xrange(track_count)]
    intercept = [x_bottom[i] - (gradient[i] * y_bottom[i]) for i in xrange(track_count)]

    # Set up track, and get pos of beginning, end points
    true_track = [ttb.Track(gradient[i], intercept[i]) for i in xrange(track_count)]
    x_true_track = [true_track[i].get_x_points() for i in xrange(track_count)]
    y_true_track = [true_track[i].get_y_points() for i in xrange(track_count)]

    # Get wire hits in each layer, then assign these to detector
    wire_hits = []
    for i in xrange(track_count): 
        wire_hits = wire_hits + (ttb.closest_hit_wires(true_detector, true_track[i]))
        
    true_detector.set_wire_hits(wire_hits)

    # New detector object for fitting
    fitting_detector = ttb.Detector(smearing=True)
    fitting_detector.set_wire_hits(wire_hits)

    # Get x-displacement of all closest approached wires
    hit_rads = [wire_hit.get_hit_dist() for wire_hit in wire_hits]
        
    wire_keys = [str((i - (i % 8)) / 8) + "-" + str(i % 8) for i in xrange(8 * track_count)] # Numbers to index all layers


    guess = [0.0 for i in xrange((2 * track_count) + 1)] # Initial guess for fitted track gradient and intercept, and module alignment (spurious convergence without this)

    print k

    fit_sigmas = [ttb.hit_resolution]*len(hit_rads)

    # Finds values for track gradient and intercept, by fitting to wire hit x-displacements
    popt, pcov = curve_fit(fitting_detector.get_hit_radius, wire_keys, hit_rads, p0=guess, sigma=fit_sigmas)
    
    alignment_residuals.append(popt[-1] - module_alignment)

    # Calculate residuals and chi-squared, to test goodness of fit.
    residuals = (np.array(hit_rads) - np.array(fitting_detector.get_hit_radius(wire_keys, popt))) / np.array(fit_sigmas)
    deg_freedom = len(wire_keys) - len(popt) # Get number of degrees of freedom for this fit (should be same for each fit)
    chi_squared = np.sum(residuals**2) #/ deg_freedom

    # Get alignment residual, and chi-squared value if reduced chi-squared less than 10
    # (Should filter out spurious fits)
    if (chi_squared < 10 * deg_freedom):
        chi_squared_vals.append(chi_squared) 
        alignment_residuals.append(popt[-1] - module_alignment)



# print ""
# print "Mean Residual:", np.mean(np.array(alignment_residuals)), "mm"
# print "Standard Deviation:", np.std(np.array(alignment_residuals)), "mm"



# # Plot histogram of residuals
# plt.hist(alignment_residuals, bins=20)
# plt.xlabel("Residual Value / mm")
# plt.ylabel("Frequency")
# plt.show()

# Plot histogram of chi-squared, then get bin contents, and edges
x_chi_plot = np.linspace(stats.chi2.ppf(0.01, deg_freedom), stats.chi2.ppf(0.99, deg_freedom), 100)
chi_bin_conts, chi_bin_edges, chi_patches = plt.hist(chi_squared_vals, bins=20, label='Calculated Chi-Squared')

# Get width of each histogram bin
chi_bin_width = (chi_bin_edges[-1] - chi_bin_edges[0]) / len(chi_bin_conts)

# Plot chi-squared distribution for given number of degrees of freedom in fit.
plt.plot(x_chi_plot, len(chi_squared_vals) * chi_bin_width * stats.chi2.pdf(x_chi_plot, deg_freedom), 'r-', label=("Chi-Squared Dist (" + str(deg_freedom) + " dof)"))
plt.legend(loc='best', frameon=False)
plt.xlabel("Chi-Squared")
plt.ylabel("Frequency")

# Create arrays of number of fits with chi-squared values within given bin, and number expected in given bin from chi-squared dist.
chi_func_bins_unfiltered = np.array([len(chi_squared_vals) * (stats.chi2.cdf(chi_bin_edges[i+1], deg_freedom) - stats.chi2.cdf(chi_bin_edges[i], deg_freedom)) for i in xrange(len(chi_bin_conts))])
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

plt.show()

