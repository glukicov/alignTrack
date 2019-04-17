# Mahalanobis 1sigma method 
# Gleb
import numpy.polynomial.polynomial as poly
import numpy as np  
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

#Define constants 
stationN=2
moduleN=8
labels=("S12", "S18")
colors=("red", "blue")

#Z-position, relative to the centre of M0, and errors [mm]
Z = [0.0, 138.4261909, 275.7030146, 412.8463143, 551.9417241, 684.9748242, 820.6609492, 956.1864672]
Z0 = Z[-1]/len(Z)
errZ = np.zeros(moduleN)

#Survey data radially [mm]
measuredR_s12 = [1.72, 2.08, 2.15, 2.29, 2.92, 2.56, 2.73, 2.92]  
measuredR_s18 = [3.05, 3.43, 3.81, 4.06, 4.46, 4.75, 5.05, 5.38]
measuredR = (measuredR_s12, measuredR_s18)
errR = 0.2 # 200 um / 0.2 mmm 
errR_fit = np.ones(moduleN)*errR # create an error array for the fit points


###Define fit function 
N=10000 # generated points 
x_fit = np.linspace(min(Z), max(Z) , N) ## arbitrarily used N=1000 points from Z_min to Z_max
# y = b*x + c 
fitEquation = r"$\displaystyle\mathrm{fit} =  b * x + c $"
def fitFunc(x, b, c):
    return b * (x) + c
    # return b * (x - Z0 ) + c


####Plot original survey data with function fit 
plt.figure(1)
plt.ylabel(r"$\Delta R$ (error="+str(errR)+") [mm]", fontsize=12)
plt.xlabel("Z-position [mm]", fontsize=12)
plt.xticks(fontsize=10, rotation=0) 
plt.yticks(fontsize=10, rotation=0)
plt.minorticks_on()
axes = plt.gca()
axes.tick_params(axis='x',which='minor',bottom=False)
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
plt.tight_layout()
plt.subplots_adjust(right=0.88)
for i_station in range(0, stationN):
    #plot survey points
    plt.plot(Z, measuredR[i_station], color=colors[i_station], label=labels[i_station], marker=".", linewidth=0)
    plt.errorbar(Z, measuredR[i_station], yerr=errR, color=colors[i_station], capsize=2, elinewidth=1, linewidth=0)
    
    ###Generate the function

    # Initial guess 
    coefs = poly.polyfit(Z, measuredR[i_station], 1, w=errR_fit, cov=True) # straight line for initial guess parameters 
    ffit = poly.polyval(x_fit, coefs) # generate y values  
    gradient, intercept=coefs[1] , coefs[0] 

    # Now fit to the custom function 
    y_func = fitFunc(x_fit, gradient, intercept)  # base function with the initial guess parameters 
    fitParams, fitCovariances = curve_fit(fitFunc, x_fit, y_func, sigma=np.ones(N)*errR)  # custom function fit with errors 
    coefs_func = (fitParams[1], fitParams[0]) # place into array (intercept, gradient)
    ffit_fucn = poly.polyval(x_fit, coefs_func)  # generate y values 

    print("gradient=", gradient, "intercept=", intercept)
    print("gradient=", fitParams[0], "intercept=", fitParams[1])

    labelStr = labels[i_station]+r": $m_z=$"+str(round(gradient,5))+r", $z_0$="+str(round(intercept,5))
    plt.plot(x_fit, ffit, color=colors[i_station])     # plot over generated points 
    plt.plot(x_fit, ffit_fucn, color="black", linestyle="--")     # plot over generated points 
    plt.text(50, 5.5-(i_station/2), labelStr)

axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
plt.savefig("Survey.png", dpi=250)