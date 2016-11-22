import toytraceback as ttb
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches
import numpy as np
import random
import sys
import math
import scipy.stats as stats
from scipy.optimize import curve_fit 

# Number of tracks to fit
track_count = 50

# Whether to include smearing of residuals, and finite straw size
detector_smear = True
detector_finite_straws = False

# Set up detector, and change alignment of first module
module_alignment = 1.5
true_detector = ttb.Detector(smearing=detector_smear, finite_straws=detector_finite_straws)
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
events = []
for j in xrange(track_count): 
    wire_hits = ttb.closest_hit_wires(true_detector, true_track[j])
    event = ttb.DetectorHitEvent(true_track[j], wire_hits, j)
    events.append(event)


# New detector object for fitting
fitting_detector = ttb.Detector(smearing=detector_smear, finite_straws=detector_finite_straws)
fitting_detector.set_events(events)

# Get x-displacement of all closest approached wires
hit_rads = []
for event in events:
    for wire_hit in event.wire_hits:
        hit_rads.append(wire_hit.get_hit_dist())

wire_keys = fitting_detector.get_hit_keys()

guess = [0.0 for i in xrange((2 * track_count) + 1)] # Initial guess for fitted track gradient and intercept, and module alignment (spurious convergence without this)

print ""
print "True Alignment:", module_alignment, "mm"
print ""
print "Fitting:"

fit_sigmas = [ttb.hit_resolution]*len(hit_rads) # Uncertainties on hit radii, from detector hit resolution

# Finds values for track gradient and intercept, by fitting to wire hit x-displacements
popt, pcov = curve_fit(fitting_detector.get_hit_radius, wire_keys, hit_rads, p0=guess, sigma=fit_sigmas)
print "Complete"
print ""

print "Fitted Alignment:", popt[-1], "mm"

# Calculate residuals and chi-squared for residuals, to test goodness of fit.
residuals = np.array(hit_rads) - np.array(fitting_detector.get_hit_radius(wire_keys, popt)) / 1.0
deg_freedom = len(residuals) - len(popt)
chi_squared = np.sum((residuals / fit_sigmas)**2)


# Print stats for all residuals
print ""
print ""
print "Hit Residual Goodness of Fit Testing:"
print ""
print "Hit Residual Mean:", np.mean(residuals), "mm"
print "Hit Residual Std-Dev:", np.std(residuals), "mm"
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
x_hit_wires = []
y_hit_wires = []


for event in fitting_detector.get_events():
    for wire_hit in event.wire_hits:
        x_hit_wires.append(wire_hit.get_wire().get_absolute_x())
        y_hit_wires.append(wire_hit.get_wire().get_absolute_y())


# Get coordinates of all wires in fitted detector
x_fitted_wires = fitting_detector.get_wires_x()
y_fitted_wires = fitting_detector.get_wires_y()


# Plot points where wires located in true and fitted detector, where wires hit in fitted detector, lines showing true track and fitted track
geom_ax = plt.subplot(1, 2, 1)
plt.scatter(x_true_wires, y_true_wires, color='red')
plt.scatter(x_fitted_wires, y_fitted_wires, color='blue')
plt.scatter(x_hit_wires, y_hit_wires, color='green', marker='*', s=100)
[plt.plot(x_true_track[i], y_true_track[i], 'r') for i in xrange(track_count)]
[plt.plot(x_fitted_track[i], y_fitted_track[i], 'g') for i in xrange(track_count)]

# Plot circles representing outside of straws, if using finite straw widths.
if (detector_finite_straws):
    for i in xrange(len(x_true_wires)):
        circle = mpatches.Circle((x_true_wires[i], y_true_wires[i]), ttb.straw_rad, facecolor='none', edgecolor='r', alpha=0.5)
        geom_ax.add_patch(circle)

# Add title
plt.title(str(track_count) + " Track Fit - Smearing: " + str(detector_smear) + " - Finite Straws: " + str(detector_finite_straws))

# Label axes, and show plot
plt.axis('equal')
plt.xlabel("x-position / mm")
plt.ylabel("y-position / mm")

hits_patch = mpatches.Patch(color='none', label=("Hits: " + str(len(residuals))))
chi_squared_patch = mpatches.Patch(color='none', label=("$\chi^2_{red}$: " + str(chi_squared / deg_freedom)))
true_align_patch = mpatches.Patch(color='none', label=("True al.: " + str(module_alignment)))
fit_align_patch = mpatches.Patch(color='none', label=("Fit al.: " + str(popt[-1])))

diag_handles = [hits_patch, chi_squared_patch, true_align_patch, fit_align_patch]
diag_labels = [hits_patch.get_label(), chi_squared_patch.get_label(), true_align_patch.get_label(), fit_align_patch.get_label()]

plt.legend(diag_handles, diag_labels, loc='best', frameon=False)

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
print "Fitted Gaussian Mean:", mu_res, "mm"
print "Fitted Gaussian Std-Dev:", sigma_res, "mm"
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

# Patches to show chi-squared, number of entries in legend
number_patch = mpatches.Patch(color='none', label=("N: " + str(sum(res_obs_bin_contents))))
chi_squared_patch = mpatches.Patch(color='none', label=("$\chi^2_{red}$: " + str(chi_squared_fitted_gauss / len(res_obs_bin_contents))))
mean_patch = mpatches.Patch(color='none', label=("$\mu$: " + str(mu_res)))
std_patch = mpatches.Patch(color='none', label=("$\sigma$: " + str(sigma_res)))

# Sort out handles and labels, then create legend, and show plot
hist_handles = [number_patch, chi_squared_patch, mean_patch, std_patch]
hist_labels = [number_patch.get_label(), chi_squared_patch.get_label(), mean_patch.get_label(), std_patch.get_label()]
plt.legend(hist_handles, hist_labels, loc='best', frameon=False)

plt.show()
