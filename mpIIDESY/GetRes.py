####################################################################
# Input: TrackerAlignment.root analysis-level TFile after alignment
# Output: Steerable by args 
#s 
# e.g. python3 ../../GetOVPUV.py -m -1 -mode plot 
# Will produce UV Residuals and Residual SD plots 
# with no (-1) modules removed 
# 
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 26 November 2018 by Gleb
#####################################################################

from ROOT import TFile, TStyle, TCanvas, gStyle, TF1
import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines 
import argparse, sys
from math import log10, floor
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
import subprocess
def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--moduleN', help='mode', default=-1) # # of removed module from tracking (if applicable)
parser.add_argument('-f', '--fileN', help='input ROOT file', default="TrackerAlignment.root")
parser.add_argument('-s', '--stationN', help='station number')
parser.add_argument('-eL', '--eL', help='extra label', default="", type=str)
args = parser.parse_args()

eL=str(args.eL)

NModules=8
NLayers=4 # per modules
NTotalLayers=32
stationN=str(args.stationN)
LayerNames = ["U0", "U1", "V0", "V1"]
moduleNamesInitial=np.arange(1, NModules+1) #1-8
layerNamesInitial=np.arange(1, NTotalLayers+1) #1-32

#Dealing with removed module 
if (int(args.moduleN) != -1):
    removedModule=int(args.moduleN)
    moduleNames=np.delete(moduleNamesInitial, removedModule-1) # indexing so -1
    removedLayers= np.arange(removedModule*4-3,removedModule*4+1)
    print("removedPlanes ", removedLayers) # Layers: 0, 1... Planes: 1, 2...
    print("layerNamesInitial", layerNamesInitial)
    layerNames=np.delete(layerNamesInitial, removedLayers-1) # indexing so -1
else:
    removedModule=-1
    moduleNames=moduleNamesInitial
    layerNames=layerNamesInitial

print("Getting Plots for", len(moduleNames), "modules: ", moduleNames, "and")
print(len(layerNames), "planes: ", layerNames)


fileName = str(args.fileN)
regime = None 
if (fileName == "TrackerAlignment.root"):
    print("Plotting Residuals from Alignment Tracks!")
    regime="align"
    # input("Correct? [press enter]")

elif (fileName == "gm2tracker_ana.root"):
    print("Plotting Residuals from Quality Geane Tracks!")
    regime="track"
    # input("Correct? [press enter]") 
else:
    print("Not expected file name!")
    sys.exit()

f = TFile.Open(fileName)
if f:
    print(str(fileName)+" is open")
else:
    print(str(fileName)+" not found")

if (regime=="align"): # only alignment data has Pz/P and implicit station number 
	label_mean = f.Get("TrackerAlignment/Hits/Labels").GetMean()
	print(label_mean)
	if(label_mean < 1280 and label_mean > 1210):
		stationN = "S12"
	if(label_mean < 1880 and label_mean > 1810):
		stationN = "S18"
	if(label_mean < 1080 and label_mean > 1010):
		stationN = "S0"
	#-------layersPz----------
	i_totalLayer=0
	yMin = 0.92
	yMax = 1.02
	plt.figure(1)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		for n in range(0, NLayers):
			i_layer=layerNames[i_totalLayer]
			name = "TrackerAlignment/UV/PzoP reduced Module " + str(i_module) + " " + str(LayerNames[n])
			#print name
			t = f.Get(str(name))
			mean = t.GetMean()
			SD = t.GetRMS()
			SDError = t.GetRMSError()
			plt.errorbar(i_layer, mean, yerr=SD, color="red") 
			plt.plot(i_layer, mean, marker="_", color="red")
			i_totalLayer+=1
	axes.set_xlim(0, i_totalLayer+1)
	axes.set_ylim(yMin, yMax)
	plt.ylabel("<Pz/P> [error = SD]")
	plt.xlabel("Layer", fontsize=10)
	plt.savefig("layersPz.png")

####### LAYERS ##############
ResidualRMS=[]
ResidualRMSError=[]
#-------LayerResiudals_Zoom----------
i_totalLayer=0
means=[]
MeanErrors=[]
yMin = -60
yMax = 60
plt.figure(71)
axes = plt.gca()
for i_module in range(0, NModules):
	line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
for i in range(0, len(moduleNames)):
    i_module=moduleNames[i]
    #for n in range(0, NLayers):
    	#i_layer=layerNames[i_totalLayer]
    if (regime=="align"):
    	name = "TrackerAlignment/Modules/Residuals UV Module " + str(i_module)
    if (regime=="track"):
    	name = "TrackSummary"+stationN+"/PerPlane/Plane"+str(i_module)+"/Measure Residuals/UVresidualsMeasPred Plane "+str(i_layer)
    t = f.Get(str(name))
    mean = t.GetMean()
    #print(i_layer, mean)
    means.append(mean*1e3)
    SD = t.GetRMS()
    SDError = t.GetRMSError()
    ResidualRMS.append(SD*1e3)
    ResidualRMSError.append(SDError*1e3)
    meanError = t.GetMeanError()
    MeanErrors.append(meanError*1e3)
    plt.errorbar(i_module, mean*1e3, yerr=meanError*1e3, color="red", markersize=15, elinewidth=3) 
    plt.plot(i_module, mean*1e3, marker="+", color="red", linewidth=3, markersize=15)
    line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
    #axes.annotate(round_sig(mean*1e3), (i_module, mean*1e3))
    i_totalLayer+=1

avgMean = sum(means)/float(len(means))
line = [[0.5,avgMean], [NModules+1, avgMean]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-")
# plt.text(32.1, avgMean, str(round_sig(avgMean)), fontsize=9)
line = [[0.5,0.0], [NModules+1, 0.0]]
#plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
axes.set_xlim(0.5, NModules+0.5)
axes.set_ylim(yMin, yMax)
plt.title("UV Residual Mean "+stationN+" "+eL, fontsize=18)
plt.ylabel("Residual Mean [um]", fontsize=18)
plt.xlabel("Module", fontsize=18)
plt.tight_layout()
plt.savefig("Residuals_L_Zoom"+str(stationN)+".png", dpi=600)

'''
#metric
U0 = np.array(means[0::4])
U1 = np.array(means[1::4])
V0 = np.array(means[2::4])
V1 = np.array(means[3::4])

dU = U0-U1
dV = V0-V1

maxdU=max(dU)
maxdV=max(dV)
meandU=np.mean(dU)
meandV=np.mean(dV)

UV_array = (U0, U1, V0, V0)

meanModule = 0
for i_module in range(0, NModules):
	for i_layer in range(0, 4):
		meanModule=meanModule+UV_array[i_layer][i_module]- ( (U0[i_module]+U1[i_module]+V0[i_module]+V1[i_module])/ 4 )
meanModule=meanModule/NTotalLayers



stringOut= str(maxdU) + " " + str(maxdV) + " " + str(meandU) + " " + str(meandV) + " " + str(np.mean(ResidualRMS)) + " " + str(meanModule)
f=open("metric_"+stationN+".txt", "w+")
f.write(stringOut)
f.close()

line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
axes.set_xlim(0.5, NTotalLayers+1)
axes.set_ylim(yMin, yMax)
plt.title("UV Residual Mean "+stationN+" "+eL, fontsize=18)
plt.ylabel("Residual Mean [um]", fontsize=18)
plt.xlabel("Layer", fontsize=18)
plt.tight_layout()
plt.savefig("Residuals_L_Zoom"+str(stationN)+".png")

#----Layer Residual SD 
i_totalLayer=0
yMin = 80
yMax = 160
means=[]
plt.figure(81)
axes = plt.gca()
for i_module in range(0, NModules):
	line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
for i in range(0, len(moduleNames)):
	i_module=moduleNames[i]
	for n in range(0, NLayers):
		i_layer=layerNames[i_totalLayer]
		plt.errorbar(i_layer, ResidualRMS[i_totalLayer], yerr=ResidualRMSError[i_totalLayer], color="red", markersize=14, elinewidth=3) 
		plt.plot(i_layer, ResidualRMS[i_totalLayer], marker="+", color="red")
		#axes.annotate(int(round_sig( ResidualRMS[i_totalLayer], 3)), (i_module,  ResidualRMS[i_totalLayer]))
		means.append(ResidualRMS[i_totalLayer])
		i_totalLayer+=1

avgMean = sum(means)/float(len(means))
line = [[0.5,avgMean], [NTotalLayers+1, avgMean]]
plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-")
#plt.text(32.1, avgMean, str(int(round_sig(avgMean))), fontsize=9)
line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
axes.set_xlim(0.5, NTotalLayers+1)
axes.set_ylim(yMin, yMax)
plt.title("UV Residual SD "+stationN+" "+eL, fontsize=20)
plt.ylabel("Residual SD /um", fontsize=18)
plt.xlabel("Layer", fontsize=20)
plt.tight_layout()
plt.savefig("ResidualsSD_L"+str(stationN)+".png")

#combine into a single .png file 
subprocess.call(["convert" , "-append", "Residuals_L_Zoom"+str(stationN)+".png" , "ResidualsSD_L"+str(stationN)+".png", "FoM_Res_"+str(stationN)+".png"])
'''
print("ROOT File analysed!")