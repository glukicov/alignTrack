import toytraceback as ttb
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
import math
import scipy.stats as stats
from scipy.optimize import curve_fit 


# Set up detector, and change alignment of first module
true_detector = ttb.Detector()
true_detector.set_module_x_align(1, 1.0)

# Get wire coordinates
x_wires = true_detector.get_wires_x()
y_wires = true_detector.get_wires_y()

# Coordinates at beginning, end of track
x_bottom = random.uniform(ttb.track_boundary_x_lower, ttb.track_boundary_x_upper)
x_top = random.uniform(ttb.track_boundary_x_lower, ttb.track_boundary_x_upper)
y_bottom = ttb.track_boundary_y_lower
y_top = ttb.track_boundary_y_upper

# Calculate track gradient and x-intercept from track coordinates
gradient = ((x_top - x_bottom) / (y_top - y_bottom))
intercept = x_bottom - (gradient * y_bottom)

# Set up track, and get pos of beginning, end points
true_track = ttb.Track(gradient, intercept)
x_true_track = true_track.get_x_points()
y_true_track = true_track.get_y_points()

print ""
print "True Params:", true_track.get_gradient(), true_track.get_intercept(), "10"
print ""

# Get wire hits in each layer, then assign these to detector 
wire_hits = ttb.closest_hit_wires(true_detector, true_track)
true_detector.set_wire_hits(wire_hits)

# New detector object for fitting
fitting_detector = ttb.Detector()
fitting_detector.set_wire_hits(wire_hits)

# Get x-displacement of all closest approached wires
x_disp = [wire_hit.get_hit_x_disp() for wire_hit in wire_hits]

layer_nums = [i for i in xrange(8)] # Numbers to index all layers

guess = [0.0, 0.0, 0.0] # Initial guess for fitted track gradient and intercept, and module alignment (spurious convergence without this)

print "Fitting:"

# Finds values for track gradient and intercept, by fitting to wire hit x-displacements
popt, pcov = curve_fit(fitting_detector.get_hit_x_displacement, layer_nums, x_disp, p0=guess)

print ""
print "Final Fitted Params:", popt
print ""
print "Final Matrix of Covariance:"
print pcov
print ""

# To draw fitted track
fitted_track = ttb.Track(popt[0], popt[1])
x_fitted_track = fitted_track.get_x_points()
y_fitted_track = fitted_track.get_y_points()

# Plot points where wires located, lines showing true track and fitted track
plt.scatter(x_wires, y_wires)
plt.plot(x_true_track, y_true_track, 'r')
plt.plot(x_fitted_track, y_fitted_track, 'g')


# Label axes, and show plot
plt.axis('equal')
plt.xlabel("x-position / mm")
plt.ylabel("y-position / mm")
plt.show()
