#
# FoM Plots for iterative runs 
#

#Header 
import argparse, sys
from scipy import stats
import matplotlib.pyplot as plt 

#Input arguments are results from runs 
parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-r', '--results', nargs='+', help='PEDE results')
args = parser.parse_args()

#Expected structure X, X_error, Y, Y_error, ...  
results=args.results
results=[float(i) for i in results]
resultsX=results[0::4]
errorsX=results[1::4]
resultsY=results[2::4]
errorsY=results[3::4]


#TODO use functions to remove redundant plotting code 
plt.figure(1)
axes=plt.gca()
run_list=[] # for int ticks on X axis 


# Radial 
plt.subplot(211) # X 
yMin = min(resultsX)*0.7 # for y-axis scaling 
yMax = max(resultsX)*1.3 
for i_run, itemX in enumerate(resultsX): 
	errorX=errorsX[i_run]
	run = i_run+1
	plt.plot(run, itemX,  marker="+", color="orange")
	plt.errorbar(run, itemX, yerr=errorX,  color="orange", markersize=12, elinewidth=1)
	run_list.append(int(run))
plt.text(1.5, yMin*1.5, "+ = Radial", fontsize=10, color="orange")
axes.set_ylim(yMin, yMax)
plt.xticks(run_list)
plt.ylabel("<|d(Recon-Truth)|>",fontsize=10)
plt.xlabel("Iteration", fontsize=10)

#Vertical 
plt.subplot(212) # Y 
yMin = min(resultsY)*0.7
yMax = max(resultsY)*1.3
for i_run, itemX in enumerate(resultsY):
	itemY=resultsY[i_run]
	errorY=errorsY[i_run]
	run = i_run+1
	plt.plot(run, itemY,  marker="+", color="purple")
	plt.errorbar(run, itemY, yerr=errorY,  color="purple", markersize=12, elinewidth=1)
plt.text(1.5, yMin*1.5, "+ : Vertical", fontsize=10, color="purple")
axes.set_ylim(yMin, yMax)
plt.xticks(run_list)
plt.ylabel("<|d(Recon-Truth)|>",fontsize=10)
plt.xlabel("Iteration", fontsize=10)

plt.subplots_adjust(hspace=.3)
plt.savefig("FoMRuns.png")


