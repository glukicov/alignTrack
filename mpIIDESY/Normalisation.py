import argparse, sys
from scipy import stats
from matplotlib.pyplot import *
import matplotlib.pyplot as plt #for plotting 
import matplotlib.ticker as ticker
import numpy as np  # smart arrays 
import itertools # smart lines
from time import gmtime, strftime 
from math import log10, floor, ceil
import subprocess, shlex 
import numpy.polynomial.polynomial as poly

#take the un-corrected truth 
T_mis_C= ( 0.0, 0.0, 0.059, 0.059, 0.004, 0.007, 0.033, -0.044, -0.064, 0.066, 0.095, 0.028, -0.071, -0.023, -0.037, -0.04)
strawModuleZPosition = [ 0.0, 134.36, 268.72, 403.08, 537.42, 671.77, 806.10, 940.406]
moduleArray=np.arange(1,9) 

for i_module in range(0, len(strawModuleZPosition)-1):
    print("dZ=", round(strawModuleZPosition[i_module]-strawModuleZPosition[i_module+1],3))

#################################
### Correct the angle in misalignment

#split into X and Y
# T_mis_C_rad = np.array(T_mis_C[0::2])
# T_mis_C_rad = [0.0, 0.059, 0.004, 0.033, -0.064, 0.095, -0.071, -0.037 ]
T_mis_C_rad = [0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0 ]
T_mis_C_ver = np.array(T_mis_C[1::2])
print("T_mis_C_rad=", T_mis_C_rad)
print("T_mis_C_ver=", T_mis_C_ver)

#form X (z-position), Y (misalignment) array 
data=[strawModuleZPosition, T_mis_C_rad] 

#Do a linear fit 
x_new = np.linspace(float(min(data[0])), float(max(data[0])), num=1e3) # generate x-points for evaluation 
coefs = poly.polyfit(data[0], data[1], 1) # straight line
ffit = poly.polyval(x_new, coefs) # plot over generated points 

print("Intercept", coefs[0])
print("Gradient", coefs[1])

#Correct the truth
corrected_T_mis_C = []
for i_module in range(0, len(moduleArray)):
    # corr = Y + Gradient * X + Intercept 
    correction = data[1][i_module] - coefs[1] * data[0][i_module] - coefs[0]
    corrected_T_mis_C.append(correction)

x_new_Corr = np.linspace(float(min(data[0])), float(max(data[0])), num=1e3) # generate x-points for evaluation 
coefs_Corr = poly.polyfit(data[0], corrected_T_mis_C, 1) # straight line
ffit_Corr = poly.polyval(x_new_Corr, coefs_Corr) # plot over generated points 

print("Corrected Intercept", coefs_Corr[0])
print("Corrected Gradient", coefs_Corr[1])
##################################
#Sanity Plot

# plt.figure(1, figsize=(12, 6)) # FoM vs Tracks
# plt.tight_layout(pad=0.4, w_pad=0.2, h_pad=1.0)
# axes = plt.subplot(111)
# plt.subplots_adjust(right=0.8)
# plt.minorticks_on()
# axes.tick_params(axis='x',which='minor',bottom=False)
# axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
# axes.set_xlim(0.9*min(strawModuleZPosition),1.1*max(strawModuleZPosition))
# axes.set_ylim(-0.140, 0.140)
# plt.plot(strawModuleZPosition, T_mis_C_rad, 'o', marker=".", color="red", label="Uncorr.")
# plt.plot(strawModuleZPosition, corrected_T_mis_C, 'o', marker="x", color="red", label="Corr.")
# axes.plot(x_new, ffit, color="red")
# axes.plot(x_new, ffit_Corr, color="red", linestyle="--")
# axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title("Misalignment X Truth Corrected")
# plt.ylabel("Misalignment [mm]")
# plt.xlabel("Module Z postion [mm]")
# plt.savefig("test.png", dpi=250)

plt.figure(1, figsize=(12, 6)) # FoM vs Tracks
plt.tight_layout(pad=0.4, w_pad=0.2, h_pad=1.0)
axes = plt.subplot(111)
plt.subplots_adjust(right=0.8)
plt.minorticks_on()
axes.tick_params(axis='x',which='minor',bottom=False)
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
axes.set_xlim(0.9*min(strawModuleZPosition),1.1*max(strawModuleZPosition))
axes.set_ylim(-0.140, 0.140)
for i_module in range(0, 8):
    rect_unc = plt.Rectangle((strawModuleZPosition[i_module]-20/2, T_mis_C_rad[i_module]-(0.05/2)), 20, 0.05, color='blue', label="Uncorr.")
    axes.add_patch(rect_unc)
    rect_cor = plt.Rectangle((strawModuleZPosition[i_module]-20/2, corrected_T_mis_C[i_module]-(0.05/2)), 20, 0.05, color='green', label="Corr.")
    axes.add_patch(rect_cor)
    rect_cor_angle = plt.Rectangle((strawModuleZPosition[i_module]-20/2, corrected_T_mis_C[i_module]-(0.05/2)), 2, 0.01, angle=30.0, color='red', label="angle.")
    axes.add_patch(rect_cor_angle)
    # plt.plot(strawModuleZPosition, T_mis_C_rad, 'o', marker=".", color="red", label="Uncorr.")
    # plt.plot(strawModuleZPosition, corrected_T_mis_C, 'o', marker="x", color="red", label="Corr.")
axes.plot(x_new, ffit, color="red")
axes.plot(x_new, ffit_Corr, color="red", linestyle="--")
axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title("Misalignment X Truth Corrected")
plt.ylabel("Misalignment [mm]")
plt.xlabel("Module Z postion [mm]")
plt.savefig("test.png", dpi=250)




# xTitle = "



# plt.plot(moduleArray, np.array(T_mis_C_rad), marker=".", color="red", label="Truth Mis. (uncorr.)")
