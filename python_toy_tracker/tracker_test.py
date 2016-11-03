import toytraceback as ttb
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
import math
import scipy.stats as stats
from scipy.optimize import curve_fit 

# Number of tracks to fit
track_count = 5

# Set up detector, and change alignment of first module
module_alignment = 0.3
true_detector = ttb.Detector()
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


print ""
print "True Params:", [str(true_track[i].get_gradient()) + " " + str(true_track[i].get_intercept()) for i in xrange(track_count)], str(module_alignment)
print ""

# Get wire hits in each layer, then assign these to detector
wire_hits = []
for j in xrange(track_count): 
    wire_hits = wire_hits + (ttb.closest_hit_wires(true_detector, true_track[j]))

true_detector.set_wire_hits(wire_hits)

# New detector object for fitting
fitting_detector = ttb.Detector()
fitting_detector.set_wire_hits(wire_hits)

# Get x-displacement of all closest approached wires
hit_rads = [wire_hit.get_hit_dist() for wire_hit in wire_hits]

wire_keys = [str((i - (i % 8)) / 8) + "-" + str(i % 8) for i in xrange(8 * track_count)] # Numbers to index all layers


guess = [0.0 for i in xrange((2 * track_count) + 1)] # Initial guess for fitted track gradient and intercept, and module alignment (spurious convergence without this)

print "Fitting:"

# Finds values for track gradient and intercept, by fitting to wire hit x-displacements
popt, pcov = curve_fit(fitting_detector.get_hit_radius, wire_keys, hit_rads, p0=guess)

print ""
print "Final Fitted Params:", popt
print ""
print "Final Matrix of Covariance:"
print pcov
print ""

# Calculate residuals and chi-squared, to test goodness of fit.
residuals = np.array(hit_rads) - np.array(fitting_detector.get_hit_radius(wire_keys, popt)) / 1.0
deg_freedom = len(wire_keys) - len(popt)
chi_squared = np.sum(residuals**2) / deg_freedom
print ""
print "Residuals:", residuals
print ""
print "Chi^2:", chi_squared
print ""

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

# Plot histogram of residuals
plt.subplot(1, 2, 2)
plt.hist(residuals, normed=True, bins=5)
plt.xlabel("Residual Value / mm")
plt.ylabel("Normalised Frequency")

plt.show()
