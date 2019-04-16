import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy.polynomial.polynomial as poly

data = [1, 2, 3, 4, 5, 6, 7, 8 ], [ 17,  2, -7, -12 ,-12 ,-7,  2 , 16]

#Do a curve x^2 fit 
x_new = np.linspace(float(min(data[0])), float(max(data[0])), num=1000) # generate x-points for evaluation 
coefs = poly.polyfit(data[0], data[1], 2) # x2 curve
ffit = poly.polyval(x_new, coefs) # plot over generated points 

print("Intercept (a)", coefs[0])
print("Gradient 1 (b)", coefs[1])
print("Gradient 2 (c)", coefs[2])

print("Midpoint", -coefs[1] / (2*coefs[2]))

data_truth = [1, 2, 3, 4, 5, 6, 7, 8 ], [12,   3 , -6 ,-10 , -9 , -5 ,  2 , 12]

#Do a curve x^2 fit 
x_new_truth = np.linspace(float(min(data_truth[0])), float(max(data_truth[0])), num=1000) # generate x-points for evaluation 
coefs_truth = poly.polyfit(data_truth[0], data_truth[1], 2) # x2 curve
ffit_truth = poly.polyval(x_new_truth, coefs_truth) # plot over generated points 

plt.figure(1) # FoM vs Tracks
plt.tight_layout(pad=0.4, w_pad=0.2, h_pad=1.0)
axes = plt.gca()
plt.subplots_adjust(right=0.65)
plt.minorticks_on()
axes.tick_params(axis='x',which='minor',bottom=False)
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
plt.plot(data[0], data[1], 'o', marker=".", color="red", label="Residuum misalignment")
plt.plot(data_truth[0], data_truth[1], 'o', marker="x", color="red", label="Truth curve")
axes.plot(x_new, ffit, color="red", linestyle="-", label="Fit: residuum misalignment")
axes.plot(x_new_truth, ffit_truth, color="red", linestyle="--", label="Fit: truth curve")
plt.ylabel("Misalignment [mm]")
plt.xlabel("Module")
axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
plt.savefig("Curve.png", dpi=250)
