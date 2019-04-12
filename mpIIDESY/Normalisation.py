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
strawModuleZPosition = [ 0.0, 134.36, 268.72, 403.08, 537.42, 671.77, 806.10, 940.406]
ModuleArray=np.arange(1,9)
ModuleN=len(ModuleArray)
globalN=2

FHICLPatchName = ["strawModuleRShift", "strawModuleHShift"] # no FHICL patch for angles in art yet  
FHICLServicePath = "services.Geometry.strawtracker." #full path of the tracking FHICL patch     print("dZ=", round(strawModuleZPosition[i_module]-strawModuleZPosition[i_module+1],3))
GlobalParNames = ["Radial", "Vertical"]
units = [r" [$\mathrm{\mu m}$]", r" [$\mathrm{\mu m}$]"]
#################################
### Correct the angle in misalignment

#split into X and Y
T_mis_C_rad = [-0.016, 0.037, -0.019, 0.014, -0.074, 0.099, -0.049, 0.008]
T_mis_C_ver = [-0.009, 0.04, -0.017, -0.067, 0.048, 0.023, -0.011, -0.006]
stationN="12"

yMax = [0.120, 0.120]
yMin = [-0.120, -0.120]

#Angles and Offsets 
corrected_truth = []
truth = [T_mis_C_rad, T_mis_C_ver]
plt.figure(1)
for i_global in range(0, globalN):
    print("Truth "+FHICLPatchName[i_global]+stationN+" :", truth[i_global])
    # Correct the truth offsets to have 0 angle and 0 overall offset via a a linear fit 
    x_new = np.linspace(float(min(strawModuleZPosition)), float(max(strawModuleZPosition)), num=1000) # generate x-points for evaluation 
    coefs = poly.polyfit(strawModuleZPosition, truth[i_global], 1) # straight line 
    ffit = poly.polyval(x_new, coefs) # plot over generated points 
    corrected_array = []
    print("Intercept", coefs[0])
    print("Gradient", coefs[1])
    for i_module in range(0, ModuleN):
        # corr = Y (in um) + Gradient * X + Intercept 
        correction =  round(truth[i_global][i_module] - coefs[1] * strawModuleZPosition[i_module]  - coefs[0],3)
        corrected_array.append(correction)
    corrected_truth.append(corrected_array)
    x_new = np.linspace(float(min(strawModuleZPosition)), float(max(strawModuleZPosition)), num=1000) # generate x-points for evaluation 
    coefs_Corr = poly.polyfit(strawModuleZPosition, corrected_truth[i_global], 1) # straight line 
    ffit_Corr = poly.polyval(x_new, coefs_Corr) # plot over generated points
    print("Corrected Intercept", coefs_Corr[0])
    print("Corrected Gradient", coefs_Corr[1])

    plt.subplot(int( str(globalN)+"1"+str(int(i_global+1)) )) 
    axes = plt.gca()
    axes.set_xlim(strawModuleZPosition[0]*0.9, strawModuleZPosition[-1]*1.1)
    axes.set_ylim(yMin[i_global], yMax[i_global])
    plt.title(GlobalParNames[i_global]+" misalignment in S"+stationN, fontsize=12)
    plt.ylabel(r"Misalignment "+ units[i_global])
    plt.xlabel("Module", fontsize=12)
    plt.xticks(fontsize=10, rotation=0) 
    plt.yticks(fontsize=10, rotation=0)
    plt.minorticks_on()
    axes.tick_params(axis='x',which='minor',bottom=False)
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
    plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
    plt.plot(strawModuleZPosition, truth[i_global], marker=".", color="red", label="Uncorr.", linewidth=0)
    plt.plot(x_new, ffit, color="red")
    plt.plot(x_new, ffit_Corr, color="red", linestyle="--")
    plt.plot(strawModuleZPosition, corrected_truth[i_global], marker="x", color="red", label="Corr.", linewidth=0)
    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 

plt.savefig("Line_Corrected.png", dpi=250)

for i_global in range(0, globalN):
     print("Corrected Truth "+FHICLPatchName[i_global]+stationN+" :", corrected_truth[i_global])

# Curve 
#Do a curve x^2 fit 

truth = [corrected_truth[0], corrected_truth[1]]
corrected_truth = []

plt.figure(2)
for i_global in range(0, globalN):
    x_new = np.linspace(float(min(strawModuleZPosition)), float(max(strawModuleZPosition)), num=1e3) # generate x-points for evaluation 
    coefs = poly.polyfit(strawModuleZPosition, truth[i_global], 2) # x2 curve
    ffit = poly.polyval(x_new, coefs) # plot over generated points 

    print("Intercept (a)", coefs[0])
    print("Gradient 1 (b)", coefs[1])
    print("Gradient 2 (c)", coefs[2])
    print("Midpoint", -coefs[1] / (2*coefs[2]))

    corrected_array = []
    for i_module in range(0, ModuleN):
        # corr = Y (in um) + Gradient * X + Intercept 
        correction =  round(truth[i_global][i_module] - coefs[2] * (strawModuleZPosition[i_module] ** 2)  - coefs[1] * strawModuleZPosition[i_module]  - coefs[0],3)
        corrected_array.append(correction)
    corrected_truth.append(corrected_array)
    print(corrected_truth)
    x_new = np.linspace(float(min(strawModuleZPosition)), float(max(strawModuleZPosition)), num=1000) # generate x-points for evaluation 
    coefs_Corr = poly.polyfit(strawModuleZPosition, corrected_truth[i_global], 1) # straight line 
    ffit_Corr = poly.polyval(x_new, coefs_Corr) # plot over generated points
    print("Corrected Intercept", coefs_Corr[0])
    print("Corrected Gradient", coefs_Corr[1])


    plt.subplot(int( str(globalN)+"1"+str(int(i_global+1)) )) 
    axes = plt.gca()
    axes.set_xlim(strawModuleZPosition[0]*0.9, strawModuleZPosition[-1]*1.1)
    axes.set_ylim(yMin[i_global], yMax[i_global])
    plt.title(GlobalParNames[i_global]+" misalignment in S"+stationN, fontsize=12)
    plt.ylabel(r"Misalignment "+ units[i_global])
    plt.xlabel("Module", fontsize=12)
    plt.xticks(fontsize=10, rotation=0) 
    plt.yticks(fontsize=10, rotation=0)
    plt.minorticks_on()
    axes.tick_params(axis='x',which='minor',bottom=False)
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
    plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
    plt.plot(strawModuleZPosition, truth[i_global], marker=".", color="red", label="Uncorr.", linewidth=0)
    plt.plot(x_new, ffit, color="red")
    plt.plot(x_new, ffit_Corr, color="red", linestyle="--")
    plt.plot(strawModuleZPosition, corrected_truth[i_global], marker="x", color="red", label="Corr.", linewidth=0)
    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 

plt.savefig("Curve_Corrected.png", dpi=250)

for i_global in range(0, globalN):
     print("Corrected Truth "+FHICLPatchName[i_global]+stationN+" :", corrected_truth[i_global])

subprocess.call(["convert" , "-append", "Line_Corrected.png" , "Curve_Corrected.png", "Corrected.png"])