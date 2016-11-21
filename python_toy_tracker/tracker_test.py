import toytraceback as ttb
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import random
import sys
import math
import scipy.stats as stats
from scipy.optimize import curve_fit 

# Number of tracks to fit
track_count = 50

# Set up detector, and change alignment of first module
module_alignment = 1.5
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


# Calculate track gradient and x-intercept from track coordinates
gradient = [((x_top[i] - x_bottom[i]) / (y_top[i] - y_bottom[i])) for i in xrange(track_count)]
intercept = [x_bottom[i] - (gradient[i] * y_bottom[i]) for i in xrange(track_count)]

# Set up track, and get pos of beginning, end points
true_track = [ttb.Track(gradient[i], intercept[i]) for i in xrange(track_count)]
x_true_track = [true_track[i].get_x_points() for i in xrange(track_count)]
y_true_track = [true_track[i].get_y_points() for i in xrange(track_count)]


# Get wire hits in each layer, then assign these to detector
wire_hits = []
for j in xrange(track_count): 
    wire_hits = wire_hits + (ttb.closest_hit_wires(true_detector, true_track[j]))

true_detector.set_wire_hits(wire_hits)

# New detector object for fitting
fitting_detector = ttb.Detector(smearing=True)
fitting_detector.set_wire_hits(wire_hits)

# Get x-displacement of all closest approached wires
hit_rads = [wire_hit.get_hit_dist() for wire_hit in wire_hits]

wire_keys = [str((i - (i % 8)) / 8) + "-" + str(i % 8) for i in xrange(8 * track_count)] # Numbers to index all layers


guess = [0.0 for i in xrange((2 * track_count) + 1)] # Initial guess for fitted track gradient and intercept, and module alignment (spurious convergence without this)

print ""
print "True Alignment:", module_alignment
print ""
print "Fitting:"

fit_sigmas = [ttb.hit_resolution]*len(hit_rads) # Uncertainties on hit radii, from detector hit resolution

# Finds values for track gradient and intercept, by fitting to wire hit x-displacements
popt, pcov = curve_fit(fitting_detector.get_hit_radius, wire_keys, hit_rads, p0=guess, sigma=fit_sigmas)
print "Complete"
print ""

print "Fitted Alignment:", popt[-1]

# Calculate residuals and chi-squared for residuals, to test goodness of fit.
residuals = np.array(hit_rads) - np.array(fitting_detector.get_hit_radius(wire_keys, popt)) / 1.0
deg_freedom = len(residuals) - len(popt)
chi_squared = np.sum((residuals / fit_sigmas)**2)


# Print stats for all residuals
print ""
print ""
print "Hit Residual Goodness of Fit Testing:"
print ""
print "Hit Residual Mean:", np.mean(residuals)
print "Hit Residual Std-Dev:", np.std(residuals)
print ""
print "Number of Hit Residuals:", len(residuals)
print "Number of Constraining Parameters", len(popt) 
print "Hit Residual Reduced chi^2:", chi_squared / deg_freedom
print "P(chi^2)", 1 - stats.chi2.cdf(chi_squared, deg_freedom)

# To draw fitted track
fitted_track = [ttb.Track(popt[2 * i], popt[(2 * i) + 1]) for i in xrange(track_count)]
x_fitted_track = [fitted_track[i].get_x_points() for i in xrange(track_count)]
y_fitted_track = [fitted_track[i].get_y_points() for i in xrange(track_count)]

# Get coordinates of wires hit in fitted detector
x_hit_wires = [wire_hit.get_wire().get_absolute_x() for wire_hit in fitting_detector.get_wire_hits()]
y_hit_wires = [wire_hit.get_wire().get_absolute_y() for wire_hit in fitting_detector.get_wire_hits()]

# Get coordinates of all wires in fitted detector
x_fitted_wires = fitting_detector.get_wires_x()
y_fitted_wires = fitting_detector.get_wires_y()


# Plot points where wires located in true and fitted detector, where wires hit in fitted detector, lines showing true track and fitted track
plt.subplot(1, 2, 1)
plt.scatter(x_true_wires, y_true_wires, color='red')
plt.scatter(x_fitted_wires, y_fitted_wires, color='blue')
plt.scatter(x_hit_wires, y_hit_wires, color='green')
[plt.plot(x_true_track[i], y_true_track[i], 'r') for i in xrange(track_count)]
[plt.plot(x_fitted_track[i], y_fitted_track[i], 'g') for i in xrange(track_count)]

# Label axes, and show plot
plt.axis('equal')
plt.xlabel("x-position / mm")
plt.ylabel("y-position / mm")

# Filter large residuals for fitting of gaussian function
filtered_residuals = [x for x in residuals if abs(x) < (10 * ttb.hit_resolution)]

# Get mean and std-dev of fitted gaussian, then array of x-values for plotting
(mu_res, sigma_res) = stats.norm.fit(filtered_residuals)
gaus_x_res = np.linspace(mu_res - 4 * sigma_res, mu_res + 4 * sigma_res, 500)


# Plot histogram of residuals
plt.subplot(1, 2, 2)
res_obs_bin_contents_unfiltered, res_obs_bin_edges, _ = plt.hist(filtered_residuals, bins=20)

# Plot fitted gaussian function
bin_width = (res_obs_bin_edges[-1] - res_obs_bin_edges[0]) / len(res_obs_bin_contents_unfiltered)
plt.plot(gaus_x_res, sum(res_obs_bin_contents_unfiltered) * bin_width * mlab.normpdf(gaus_x_res, mu_res, sigma_res), 'r')
plt.xlabel("Residual Value / mm")
plt.ylabel("Frequency")

# Output Gaussian parameters
print ""
print "Gaussian Goodness-of-Fit Testing:"
print ""
print "Fitted Gaussian Mean:", mu_res
print "Fitted Gaussian Std-Dev:", sigma_res
print ""

# Get expected number of residuals in each bin, from fitted gaussian
bin_count = len(res_obs_bin_contents_unfiltered)
res_func_bin_contents_unfiltered = [len(filtered_residuals) * (stats.norm.cdf(res_obs_bin_edges[i+1], scale=(sigma_res), loc=mu_res) - stats.norm.cdf(res_obs_bin_edges[i], scale=(sigma_res), loc=mu_res)) for i in xrange(bin_count)]

# Filter out bins containing zero observed residuals, to exclude from chi-squared calculation
res_fun_bin_contents = np.array([res_func_bin_contents_unfiltered[i] for i in xrange(bin_count) if res_obs_bin_contents_unfiltered[i] > 0])
res_obs_bin_contents = np.array([res_obs_bin_contents_unfiltered[i] for i in xrange(bin_count) if res_obs_bin_contents_unfiltered[i] > 0])

res_sigmas = np.sqrt(res_fun_bin_contents) # Uncertainty in observed bins (Poisson)

# Calculate chi-squared for histogram of residuals and gaussian fit
chi_squared_fitted_gauss = np.sum(((res_fun_bin_contents - res_obs_bin_contents) / res_sigmas)**2)

print "Number of bins used for test:", len(res_obs_bin_contents)
print "Reduced Chi-squared of measured residuals to fitted gaussian:", chi_squared_fitted_gauss / len(res_obs_bin_contents)
print "P(chi^2):", 1 - stats.chi2.cdf(chi_squared_fitted_gauss, len(res_obs_bin_contents))
print ""

plt.show()
