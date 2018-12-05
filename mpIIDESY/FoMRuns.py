####################################################################
# FoM Plots for iterative runs  
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 26 November 2018 by Gleb
#####################################################################

#Header 
import argparse, sys
from scipy import stats
import matplotlib.pyplot as plt 
from matplotlib.ticker import FormatStrFormatter # ensure floats as tick labels
import matplotlib.ticker as ticker
import numpy as np

#Input arguments are results from runs 
parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-r', '--results', nargs='+', help='PEDE results')
args = parser.parse_args()

#Expected structure X, X_error, Y, Y_error, ...  
results=args.results
results=[float(i) for i in results]
# resultsX=results[0::4]
# errorsX=results[1::4]
# resultsY=results[2::4]
# errorsY=results[3::4]

resultsY=results[0::2]
errorsY=results[1::2]

resultsY=np.diff(resultsY)

print(resultsY)
print(errorsY)

#TODO use functions to remove redundant plotting code 
plt.figure(1)
axes=plt.gca()
run_list=[] # for int ticks on X axis 


# # Radial 
# plt.subplot(211) # X 
# yMin = min(resultsX)*0.7 # for y-axis scaling 
# yMax = max(resultsX)*1.3 
# for i_run, itemX in enumerate(resultsX): 
# 	errorX=errorsX[i_run]
# 	run = i_run+1
# 	plt.plot(run, itemX,  marker="+", color="orange")
# 	plt.errorbar(run, itemX, yerr=errorX,  color="orange", markersize=12, elinewidth=1)
# 	run_list.append(int(run))
# plt.text(1.5, yMin*1.5, "+ = Radial", fontsize=10, color="orange")
# axes.set_ylim(yMin, yMax)
# plt.xticks(run_list)
# plt.ylabel("<|d(Recon-Truth)|>",fontsize=10)
# plt.xlabel("Iteration", fontsize=10)

# #Vertical 
# plt.subplot(212) # Y 
# yMin = min(resultsY)*0.7
# yMax = max(resultsY)*1.3
# for i_run, itemX in enumerate(resultsY):
# 	itemY=resultsY[i_run]
# 	errorY=errorsY[i_run]
# 	run = i_run+1
# 	plt.plot(run, itemY,  marker="+", color="purple")
# 	plt.errorbar(run, itemY, yerr=errorY,  color="purple", markersize=12, elinewidth=1)
# plt.text(1.5, yMin*1.5, "+ : Vertical", fontsize=10, color="purple")
# axes.set_ylim(yMin, yMax)
# plt.xticks(run_list)
# plt.ylabel("<|d(Recon-Truth)|>",fontsize=10)
# plt.xlabel("Iteration", fontsize=10)

# plt.subplots_adjust(hspace=.3)
# plt.savefig("FoMRuns.png")


#Mis 
yMin = min(resultsY)*0.4
yMax = max(resultsY)*1.6
for i_run, itemX in enumerate(resultsY):
	itemY=resultsY[i_run]
	errorY=errorsY[i_run]
	run = i_run+2
	run_list.append(int(run))
	plt.plot(run, itemY,  marker="+", color="purple")
	plt.errorbar(run, itemY, yerr=errorY,  color="purple", markersize=12, elinewidth=2)
#plt.text(1.5, yMin*1.5, "+ : Vertical", fontsize=10, color="purple")
axes.set_ylim(yMin, yMax)
plt.title("Discrete difference vs iteration", fontsize=16)
plt.ylabel("<|d(ReconX)|>",fontsize=16)
plt.xlabel("Iteration", fontsize=16)
# axes.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True)) # ensure floats as tick labels
# plt.subplots_adjust(hspace=.3)
plt.xticks(run_list)
plt.savefig("FoMRuns.png")



