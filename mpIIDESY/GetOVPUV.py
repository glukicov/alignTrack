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
parser.add_argument('-m', '--moduleN', help='mode', default=-1)
parser.add_argument('-f', '--fileN', help='input ROOT file', default="TrackerAlignment.root")
parser.add_argument("-mode", "--mode", default="plot")
parser.add_argument("-pvals", "--pvals", nargs='+')
parser.add_argument("-p0", "--p0", nargs='+')
parser.add_argument("-mis", "--mis", nargs='+')
parser.add_argument("-abs", "--abs")
parser.add_argument('-s', '--stationN', help='station number')
args = parser.parse_args()


NModules=8
NLayers=4 # per modules
NTotalLayers=32
mode = str(args.mode)
stationN = "S"+str(args.stationN)
LayerNames = ["U0", "U1", "V0", "V1"]
moduleNamesInitial=np.arange(1,NModules+1) #1-8
layerNamesInitial=np.arange(1, NTotalLayers+1) #1-32


#Metric

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

if (mode == "plot"):

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

	if (regime=="align"): # only alignment data has Pz/P     
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
	#
	i_totalLayer=0
	means=[]
	MeanErrors=[]
	yMin = -65
	yMax = 80
	plt.figure(71)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		for n in range(0, NLayers):
			i_layer=layerNames[i_totalLayer]
			if (regime=="align"):
				name = "TrackerAlignment/UV/Residuals UV Module " + str(i_module) + " " + str(LayerNames[n])
			if (regime=="track"):
				name = "TrackSummary"+stationN+"/PerPlane/Plane"+str(i_layer)+"/Measure Residuals/UVresidualsMeasPred Plane "+str(i_layer)
			t = f.Get(str(name))
			mean = t.GetMean()
			means.append(mean*1e3)
			SD = t.GetRMS()
			SDError = t.GetRMSError()
			ResidualRMS.append(SD*1e3)
			ResidualRMSError.append(SDError*1e3)
			meanError = t.GetMeanError()
			MeanErrors.append(meanError*1e3)
			plt.errorbar(i_layer, mean*1e3, yerr=meanError*1e3, color="red", markersize=14, elinewidth=3) 
			plt.plot(i_layer, mean*1e3, marker="+", color="red")
			#axes.annotate(round_sig(mean*1e3), (i_module, mean*1e3))
			i_totalLayer+=1

	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NTotalLayers+1, avgMean]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-")
	plt.text(32.1, avgMean, str(round_sig(avgMean)), fontsize=9)
	for i in range(0, len(means)):
		number = (means[i]-avgMean)/(MeanErrors[i])
		if (number != 0):
			number = number
		else:
			number = 0.0
		#axes.annotate( "("+str(round_sig(number,2))+")", (i+1-0.4, -0.014*1e3))

	line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Residual Mean "+stationN, fontsize=18)
	plt.ylabel("Residual Mean [um]", fontsize=18)
	plt.xlabel("Layer", fontsize=18)
	plt.tight_layout()
	plt.savefig("Residuals_L_Zoom.png")

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
	plt.text(32.1, avgMean, str(int(round_sig(avgMean))), fontsize=9)
	line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Residual SD "+stationN, fontsize=20)
	plt.ylabel("Residual SD /um", fontsize=18)
	plt.xlabel("Layer", fontsize=20)
	plt.tight_layout()
	plt.savefig("ResidualsSD_L.png")

	#combine into a single .png file 
	subprocess.call(["convert" , "-append", "Residuals_L_Zoom.png" , "ResidualsSD_L.png", "FoM_Res_"+str(stationN)+".png"])
	
	print("ROOT File analysed!")

pvals=[]
#-----pVal vs iteration----
#
if (mode == "pVal"):
	pvals=args.pvals
	pVals=[]
	#errors=[]
	pVals=args.pvals
	#errors=args.pvals[1::2]
	print(pVals)
	# for i in range (0, trialN):
	# 	pVals.append(args.pvals[i*2])	
	# 	errors.append(args.pvals[i*2+1])
	trialN = len(pVals)

	yMin = 0.001
	yMax = 0.01
	plt.figure(1)
	axes = plt.gca()

	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')

	for i in range(0, trialN):
		#plt.errorbar(int(i), float(pVals[i]), yerr=float(errors[i]), color="red") 
		plt.plot(int(i), float(pVals[i]), marker="*", color="red")

	axes.set_xlim(-0.5, trialN-0.5)
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.set_ylim(0.2, 0.6)
	line = [[-0.5,0.5], [8+0.5, 0.5]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'purple')
	#plt.ylabel(r'$\chi^{2}/ndf$ of the fit to p-value dist.', fontsize=20)
	plt.ylabel("Mean p-value", fontsize=20)
	#plt.xlabel("Module Removed", fontsize=20)
	plt.xlabel("Iteration", fontsize=20)
	plt.savefig("pFoM.png")

	print("pVal FoM produced!")

p0=[]
#-----p0 vs iteration----
#
if (mode == "p0"):
	p0=args.p0
	print(p0)
	# for i in range (0, trialN):
	# 	pVals.append(args.pvals[i*2])	
	# 	errors.append(args.pvals[i*2+1])
	trialN = len(p0)

	yMin = 0.003
	yMax = 0.01
	plt.figure(1)
	axes = plt.gca()

	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')

	for i in range(0, trialN):
		#plt.errorbar(int(i), float(pVals[i]), yerr=float(errors[i]), color="red") 
		plt.plot(int(i), float(p0[i]), marker="*", color="red")

	axes.set_xlim(-0.5, trialN-0.5)
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.set_ylim(yMin, yMax)
	plt.ylabel("p0 of the fit "+stationN, fontsize=20)
	plt.xlabel("Module Removed", fontsize=20)
	# plt.xlabel("Iteration", fontsize=20)
	plt.savefig("p0FoM.png")

	print("p0 FoM produced!")

#-----Truth Misalignment----
#
if (mode == "mis"):
	misXStr=args.mis[0:8]
	misYStr=args.mis[8:16]
	print("X Mis:: ", misXStr)
	print("Y Mis:: ",misYStr)

	#absolute misalignment
	if (args.abs == "Y"):
		for i in range(0, len(misXStr)):
			misXStr[i] = abs(int(round(float(misXStr[i])*1e3)))
			misYStr[i] = abs(int(round(float(misYStr[i])*1e3)))
	else:
		#str in mm to float in um 
		for i in range(0, len(misXStr)):
			misXStr[i] = int(round(float(misXStr[i])*1e3))
			misYStr[i] = int(round(float(misYStr[i])*1e3))


	misX = np.array(misXStr)
	misY = np.array(misYStr)

	yMin = -400
	yMax = 400
	plt.subplot(211)
	axes = plt.gca()
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		plt.plot(i_module, float(misX[i_module-1]), marker="+", color="red")
		#mod = read_png('mod.png')
		
		#Add tracker image 
		# imagebox = OffsetImage(mod, zoom=0.15)
		# xy = (i_module-0.5, float(misX[i_module-1])+150) # coordinates to position this image
		# ab = AnnotationBbox(imagebox, xy, xybox=(30., -30.),  xycoords='data', boxcoords="offset points", frameon=False)                                  
		# axes.add_artist(ab)
		axes.annotate(misX[i_module-1], (i_module-0.3, float(misX[i_module-1])+20))
	
	avgMean = sum(misX)/float(len(misX))
	SD = np.std(np.array(misX))
	line = [[0.5,avgMean], [NModules+1.5, avgMean]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'black', linestyle="-")
	plt.text(9.1, avgMean, str(int(round(avgMean))), fontsize=9)

	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+1)
	axes.set_ylim(yMin, yMax)
	plt.title("Misalignment X (<X> = %s um $\sigma$= %s um)"+stationN %(int(round(avgMean)), int(round(SD))), fontsize=14)
	plt.ylabel("Misalignment [um]", fontsize=16)

	plt.subplot(212)
	axes2 = plt.gca()
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		plt.plot(i_module, float(misY[i_module-1]), marker="+", color="red")
		#mod = read_png('mod.png')
		
		#Add tracker image 
		# imagebox = OffsetImage(mod, zoom=0.15)
		# xy = (i_module-0.5, float(misY[i_module-1])+120) # coordinates to position this image
		# ab = AnnotationBbox(imagebox, xy, xybox=(30., -30.),  xycoords='data', boxcoords="offset points", frameon=False)                                  
		# axes2.add_artist(ab)
		axes2.annotate(misY[i_module-1], (i_module-0.3, float(misY[i_module-1])+20))

	avgMean = sum(misY)/float(len(misY))
	SD = np.std(np.array(misY))
	line = [[0.5,avgMean], [NModules+1.5, avgMean]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'black', linestyle="-")
	plt.text(9.1, avgMean, str(int(round(avgMean))), fontsize=9)
	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes2.set_xlim(0.5, NModules+1)
	axes2.set_ylim(yMin, yMax)
	plt.title("Misalignment Y (<Y> = %s um $\sigma$= %s um)"+stationN %(int(round(avgMean)), int(round(SD))), fontsize=14)
	plt.ylabel("Misalignment [um]", fontsize=16)
	plt.xlabel("Module", fontsize=16)

	#plt.show()
	plt.savefig("Misalignment.png")

	print("Misalignment FoM produced!")


