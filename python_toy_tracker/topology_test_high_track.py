import toytraceback as ttb
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches
import numpy as np
import random
import sys
import math
import scipy.stats as stats
import time
from scipy.optimize import curve_fit 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import LogNorm, Normalize

start_time = time.time()

# Number of tracks to fit
track_count = 5

# Whether to include smearing of residuals, and finite straw size
detector_smear = True
detector_finite_straws = True

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

fit_sigmas = [ttb.hit_resolution]*len(hit_rads) # Uncertainties on hit radii, from detector hit resolution

popt, pcov = curve_fit(fitting_detector.get_hit_radius, wire_keys, hit_rads, p0=guess, sigma=fit_sigmas, method='trf')

# Calculate residuals and chi-squared for residuals, to test goodness of fit.
residuals = np.array(hit_rads) - np.array(fitting_detector.get_hit_radius(wire_keys, popt)) / 1.0
deg_freedom = len(residuals) - len(popt)
chi_squared = np.sum((residuals / fit_sigmas)**2)

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

fig = plt.figure()
ax = fig.add_subplot(111)

# Plot points where wires located in true and fitted detector, where wires hit in fitted detector, lines showing true track and fitted track
plt.scatter(x_true_wires, y_true_wires, color='red')
plt.scatter(x_fitted_wires, y_fitted_wires, color='blue')
plt.scatter(x_hit_wires, y_hit_wires, color='green', marker='*', s=100)
[plt.plot(x_true_track[i], y_true_track[i], 'r') for i in xrange(track_count)]
[plt.plot(x_fitted_track[i], y_fitted_track[i], 'g') for i in xrange(track_count)]

# Plot circles representing outside of straws, if using finite straw widths.
if (detector_finite_straws):
    for i in xrange(len(x_true_wires)):
        circle = mpatches.Circle((x_true_wires[i], y_true_wires[i]), ttb.straw_rad, facecolor='none', edgecolor='r', alpha=0.5)
        ax.add_patch(circle)

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
plt.show()
plt.clf()


#############################################
# True Params
#############################################
# 1D Topology
#############################################

# Test variation of chi-squared with changing alignment values, at true and 
# fitted values of track gradient and intercept

# Number of alignment values to test
test_count_1d = 500

# Array to contain calculated chi-squared values at each alignment value
chi_sq_true = []
chi_sq_fit = []

# Array to contain alignment values to be tested, and separation between test values
test_aligns = []
test_separation = ((3.0) - (-3.0)) / test_count_1d

# Create array of parameters for true values of track gradients, intercepts, and module alignment
params_true = []
for j in xrange(len(gradient)):
    params_true.append(gradient[j])
    params_true.append(intercept[j])
params_true.append(module_alignment)


# Loop through test values
for i in xrange(test_count_1d):

    
    # Get test alignment value, then get hit radii with this value
    test_align = -3.0 + (i * test_separation)
    test_aligns.append(test_align)

    # Set last entry in parameters to test alignment value. Calculate hit radii
    params_true[-1] = test_align
    test_hit_rads_true = fitting_detector.get_hit_radius(wire_keys, params_true)

    # Calculate chi-squared value, and append to array
    chi_sq_true.append(np.sum(((np.array(hit_rads) - np.array(test_hit_rads_true)) / np.array(fit_sigmas))**2) / len(wire_keys))

    # Get fitted parameters, and set final parameter to test alignment. Calculate hit radii
    params_fit = popt.copy()
    params_fit[-1] = test_align
    test_hit_rads_fit = fitting_detector.get_hit_radius(wire_keys, params_fit)

    # Calculate chi-squared value, and append to array
    chi_sq_fit.append(np.sum(((np.array(hit_rads) - np.array(test_hit_rads_fit)) / np.array(fit_sigmas))**2) / len(wire_keys))


# Plot graph chi-squared values with log scale
plt.plot(test_aligns, chi_sq_true, color='blue', label='True Track Params')
plt.plot(test_aligns, chi_sq_fit, color='red', label='Fitted Track Params')
plt.yscale('log')
plt.legend(loc='best')
plt.grid()

plt.xlabel("Module Alignment")
plt.ylabel("$\chi^2_{red}$")
plt.title("Variation of $\chi^2_{red}$ With Module Alignment, \n for True and Fitted Track Parameters")

plt.show()
plt.clf()


##############################################
# 2D Topology
##############################################

# Test variation of chi-squared with changing track gradient and intercept
# values, at true and fitted values of module alignment, and other track params

# Number of values of gradient and intercept to test
test_count_2d = 100

# 2D arrays to represent chi-squared surfaces with varying gradient and intercept
# For arbitrary track, with true and fitted parameters
chi_sq_surf_true = np.zeros((test_count_2d, test_count_2d))
chi_sq_surf_fit = np.zeros((test_count_2d, test_count_2d))

# Calculate maximum, minimum possible gradients
max_grad = (ttb.track_boundary_x_upper - ttb.track_boundary_x_lower) / (ttb.track_boundary_y_upper - ttb.track_boundary_y_lower)
min_grad = (ttb.track_boundary_x_lower - ttb.track_boundary_x_upper) / (ttb.track_boundary_y_upper - ttb.track_boundary_y_lower)

# Calculate separations between test gradients and intercepts
grad_separation = (max_grad - min_grad) / test_count_2d
int_separation = (ttb.track_boundary_x_upper - ttb.track_boundary_x_lower) / test_count_2d

# Generate arrays of test gradient, intercept values, then turn into mesh arrays
grad_array = np.arange(min_grad, max_grad, grad_separation)
int_array = np.arange(ttb.track_boundary_x_lower, ttb.track_boundary_x_upper, int_separation)
grad_array, int_array = np.meshgrid(grad_array, int_array)

# Generate array of true parameters
params_true = []
for i in xrange(len(gradient)):
    params_true.append(gradient[i])
    params_true.append(intercept[i])
params_true.append(module_alignment)

# Get array of fit parameters
params_fit = popt.copy()
fit_grad = popt[0]
fit_int = popt[1]

# Loop across test values for gradient, intercept
for i in xrange(test_count_2d):

    print i

    for j in xrange(test_count_2d):

        # Set entries in arrays of parameters to test values
        params_true[0] = grad_array[i][j]
        params_true[1] = int_array[i][j]
        params_fit[0] = grad_array[i][j]
        params_fit[1] = int_array[i][j]

        # Get arrays of hit radii
        test_hit_rads_true = fitting_detector.get_hit_radius(wire_keys, params_true)
        test_hit_rads_fit = fitting_detector.get_hit_radius(wire_keys, params_fit)
        
        # Calculate value of reduced chi-squared at this point, and add to array
        chi_sq_surf_true[i][j] = np.sum(((np.array(hit_rads) - np.array(test_hit_rads_true)) / np.array(fit_sigmas))**2 / len(wire_keys))
        chi_sq_surf_fit[i][j] = np.sum(((np.array(hit_rads) - np.array(test_hit_rads_fit)) / np.array(fit_sigmas))**2 / len(wire_keys))
 

minloc_true = (grad_array[np.unravel_index(chi_sq_surf_true.argmin(), chi_sq_surf_true.shape)], int_array[np.unravel_index(chi_sq_surf_true.argmin(), chi_sq_surf_true.shape)]) 
minloc_fit = (grad_array[np.unravel_index(chi_sq_surf_fit.argmin(), chi_sq_surf_fit.shape)], int_array[np.unravel_index(chi_sq_surf_fit.argmin(), chi_sq_surf_fit.shape)]) 


grad_array = grad_array - (grad_separation / 2)
int_array = int_array - (int_separation / 2)


#       
# True values
#
# Set up figure, axis, then get colour map
cmap = plt.get_cmap('inferno')

# Draw surface as 2D image, then show colour bar
im = plt.pcolormesh(grad_array, int_array, chi_sq_surf_true, norm=LogNorm(vmin=chi_sq_surf_true.min(), vmax=chi_sq_surf_true.max()), cmap=cmap)
plt.colorbar(im)


# Show point in parameter-space where surface minimum, and true values are
plt.plot(gradient[0], intercept[0], 'ro', label='True Values')
plt.plot(minloc_true[0], minloc_true[1], 'bo', label='Surface Minimum')

plt.legend(loc='best')
plt.xlabel("Track Gradient")
plt.ylabel("Track Intercept")
plt.title("Effect on $\chi^2_{red}$ of Varying Parameters for Single Track, \n with True Values for Other Parameters")

# Show plot, then clear
plt.show()
plt.clf()

#
# Fitted Values
#
# Draw surface as 2D image, then show colour bar
im = plt.pcolormesh(grad_array, int_array, chi_sq_surf_fit, norm=LogNorm(vmin=chi_sq_surf_fit.min(), vmax=chi_sq_surf_fit.max()), cmap=cmap)
plt.colorbar(im)


# Show point in parameter-space where surface minimum, and fitted values are
plt.plot(gradient[0], intercept[0], 'ro', label='Fitted Values')
plt.plot(minloc_fit[0], minloc_fit[1], 'bo', label='Surface Minimum')

plt.legend(loc='best')
plt.xlabel("Track Gradient")
plt.ylabel("Track Intercept")
plt.title("Effect on $\chi^2_{red}$ of Varying Parameters for Single Track, \n with Fitted Values for Other Parameters")

# Show plot, then clear
plt.show()
plt.clf()



