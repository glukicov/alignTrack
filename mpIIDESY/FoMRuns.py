import argparse, sys
from scipy import stats
import matplotlib.pyplot as plt 

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-r', '--results', nargs='+', help='PEDE results')
args = parser.parse_args()


results=args.results
results=[float(i) for i in results]
resultsX=results[0::4]
errorsX=results[1::4]
resultsY=results[2::4]
errorsY=results[3::4]

yMin = 0
yMax = 150
plt.figure(1)
axes=plt.gca()
run_list=[]
for i_run, itemX in enumerate(resultsX):
	itemY=resultsY[i_run]
	errorX=errorsX[i_run]
	errorY=errorsY[i_run]
	run = i_run+1

	plt.plot(run, itemX,  marker="+", color="orange")
	plt.errorbar(run, itemX, yerr=errorX,  color="orange", markersize=12, elinewidth=1)
	plt.plot(run, itemY,  marker="+", color="purple")
	plt.errorbar(run, itemY, yerr=errorY,  color="purple", markersize=12, elinewidth=1)

	run_list.append(int(run))


plt.text(1.5, 140, "Radial", fontsize=8, color="orange")
plt.text(1.5, 130, "Vertical", fontsize=8, color="purple")


axes.set_ylim(yMin, yMax)
plt.xticks(run_list)
plt.ylabel("<|d(Recon-Truth)|>",fontsize=8)
plt.xlabel("Iteration", fontsize=8)
plt.savefig("FoMRuns.png",  dpi=300)


