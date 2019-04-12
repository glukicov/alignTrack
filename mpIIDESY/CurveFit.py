import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy.polynomial.polynomial as poly

data = [1, 2, 3, 4, 5, 6, 7, 8 ], [ -0.032, 0.034, -0.011, 0.026, -0.062, 0.105, -0.051, -0.008]

#Do a curve x^2 fit 
x_new = np.linspace(float(min(data[0])), float(max(data[0])), num=1e3) # generate x-points for evaluation 
coefs = poly.polyfit(data[0], data[1], 2) # x2 curve
ffit = poly.polyval(x_new, coefs) # plot over generated points 

print("Intercept (a)", coefs[0])
print("Gradient 1 (b)", coefs[1])
print("Gradient 2 (c)", coefs[2])

print("Midpoint", -coefs[1] / (2*coefs[2]))

plt.figure(1, figsize=(12, 6)) # FoM vs Tracks
plt.tight_layout(pad=0.4, w_pad=0.2, h_pad=1.0)
axes = plt.subplot(111)
plt.subplots_adjust(right=0.8)
plt.minorticks_on()
axes.tick_params(axis='x',which='minor',bottom=False)
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
plt.plot(data[0], data[1], 'o', marker=".", color="red", label="Uncorr.")
axes.plot(x_new, ffit, color="blue", linestyle="--")
plt.title("Curve fit")
plt.ylabel("Misalignment [mm]")
plt.xlabel("Module #")
plt.savefig("Curve.png", dpi=250)