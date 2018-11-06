#!/usr/bin/python

####################################################################
# Sanity plots for Tracker Alignment.  
# FoM Plots for comparison of actual misalignment vs PEDE results 
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 8 May 2018 by Gleb
#####################################################################
import argparse, sys
import ROOT as r
from ROOT import *
from scipy import stats

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--mode', help='mode')
parser.add_argument('-eL', '--eL', help='label')
parser.add_argument('-s', '--stationN', help='station number')

args = parser.parse_args()
from math import log10, floor, ceil

def round_sig(x, sig=2):
	return round(x, sig-int(floor(log10(abs(x))))-1)

def float_round(num, places = 0, direction = floor):
    return direction(num * (10**places)) / float(10**places)

extraLabel =-1
extraLabel=str(args.eL)
stationN = str(args.stationN)

file = str(args.mode)
from matplotlib.pyplot import *
import matplotlib.pyplot as plt #for plotting 
import matplotlib.ticker as ticker
import numpy as np  # smart arrays 
import itertools # smart lines
from time import gmtime, strftime 
import subprocess
import plotly.plotly
import plotly.graph_objs as go
import pandas as pd
plotly.tools.set_credentials_file(username='glebluk', api_key='FK1MEM1aDROhONaqC7v7')


#Truth Misalignment 

# T_mis_C=(0.1, -0.07, -0.08, 0.05, 0.15, 0.1, -0.04, 0.01) # Case C (Initial)
# expectPars = (11, 21, 31, 41, 51, 61, 71, 81)
if (stationN == "10"):
	expectPars = (1011, 1012, 1021, 1022, 1031, 1032, 1041, 1042, 1051, 1052, 1061, 1062, 1071, 1072, 1081, 1082)  # XY  [Station 0]

if (stationN == "18"):
	expectPars = (1811, 1812, 1821, 1822, 1831, 1832, 1841, 1842, 1851, 1852, 1861, 1862, 1871, 1872, 1881, 1882)  # XY  [Station 18]

if (stationN == "12"):
	expectPars = (1211, 1212, 1221, 1222, 1231, 1232, 1241, 1242, 1251, 1252, 1261, 1262, 1271, 1272, 1281, 1282)  # XY  [Station 12]

#Truth Misalignment 



# T_mis_C=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) # No Misalignment

# MUSE 
if (stationN == "12"):
	T_mis_C = (0.0, 0.0, 0.0, 0.0, -0.2, 0.0, 0.1, 0.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ) #S12 MUSE 
if (stationN == "18"):
	T_mis_C = (0.0, 0.0, 0.0, 0.0, -0.2, -0.1, 0.1, 0.0, 0.15, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0 ) #S18 MUSE 

# #DATA 
# if (stationN == "12"):
# 	T_mis_C = (0.063, -0.012, -0.125, 0.015, -0.219, 0.091, -0.207, 0.053, -0.22, 0.013, -0.152, -0.03, -0.029, -0.041, 0.09, 0.044 ) #S12 16367 
# if (stationN == "18"):
# 	T_mis_C = (-0.069, 0.013, -0.02, -0.006, 0.052, 0.157, 0.143, 0.154, 0.08, 0.156, 0.099, 0.212, 0.095, 0.094, -0.052, -0.115) #S18 16367 


globalN=int(len(expectPars))/int(8)
offsets = [0 for i in xrange(8*globalN)] # first set to 0 unless doing iterations
useOffsets = False 

# Run 0 
mis_C=T_mis_C  # the truth is the only misalignment 
print "Initial Truth Misalignment [mm]: ", mis_C
print "With expected Parameters: ", expectPars
raw_input("Truth Misalignment correct? [press enter]") 

# ----------------------------
#Run 1
useOffsets = True
##MUSE 
# if (stationN == "12"):
# 	offsets = ( 0.03, 0.009, -0.007, -0.007, -0.158, -0.002, 0.026, -0.006, 0.055, -0.006, -0.033, -0.001, -0.014, -0.003, 0.001, 0.004 ) # S12 Run 1
# if (stationN == "18"):
# 	offsets = ( 0.031, 0.005, -0.007, -0.004, -0.159, -0.055, 0.025, -0.004, 0.055, -0.008, -0.033, -0.036, -0.014, -0.003, 0.001, 0.004 ) # S18 Run 1

# if (stationN == "12"):
# 	offsets = (  -0.02, -0.003, 0.018, 0.003, -0.093, -0.002, 0.025, -0.005, 0.03, 0.003, -0.029, -0.0, -0.001, 0.001, 0.001, -0.002 ) # S12 Run 1
# if (stationN == "18"):
# 	offsets = (  -0.02, -0.017, 0.017, 0.015, -0.093, -0.047, 0.025, -0.002, 0.029, 0.005, -0.029, -0.022, -0.001, 0.003, 0.001, -0.004 ) # S18 Run 1

# offsets = (0.053, 0.01, -0.026, -0.009, -0.216, -0.051, -0.008, 0.008, 0.025, -0.001, -0.071, -0.033, -0.027, 0.002, 0.018, -0.002) # S18 Run 1
# offsets = (0.055, 0.014, -0.01, -0.012, -0.238, -0.092, 0.028, -0.011, 0.072, -0.014, -0.059, -0.049, -0.028, -0.005, 0.006, 0.007) # S18 Run 2
# offsets = (0.062, -0.001, 0.003, 0.001, -0.233, -0.074, 0.045, -0.0, 0.09, -0.003, -0.051, -0.034, -0.029, -0.001, -0.007, 0.003) # S18 Run 3
# offsets = (0.065, 0.003, 0.019, -0.002, -0.213, -0.081, 0.068, 0.005, 0.112, 0.002, -0.036, -0.035, -0.027, 0.001, -0.022, 0.001) # S18 Run 4
offsets = (0.072, 0.005, 0.032, -0.004, -0.196, -0.082, 0.089, -0.004, 0.13, -0.004, -0.025, -0.037, -0.027, -0.003, -0.037, 0.006) # S18 Run 5

print "Offsets Run 1 [mm]: ", offsets
raw_input("Offsets :: Run 1 correct? [press enter]") 
# ----------------------------


# # # ----------------------------
# #Run 2
# useOffsets = True
# ##MUSE 
# if (stationN == "12"):
# 	offsets = (  0.031, -0.003, 0.012, 0.004, -0.179, -0.006, 0.081, 0.002, 0.121, -0.003, -0.012, 0.006, -0.01, 0.004, -0.019, -0.003 ) # S12 Run 1
# if (stationN == "18"):
# 	offsets = (  0.023, 0.009, 0.018, -0.007, -0.16, -0.058, 0.108, -0.003, 0.145, -0.004, 0.008, -0.045, -0.003, -0.004, -0.027, 0.006 ) # S18 Run 1

# print "Offsets Run 2 [mm]: ", offsets
# raw_input("Offsets :: Run 2 correct? [press enter]") 
# # # ----------------------------

if ( len(mis_C) != len(expectPars) ):
	print "Enter Truth data in the right format!"

print "Truth Misalignments and offsets [mm]: "
print ["{0:0.3f}".format(i) for i in mis_C]

print "Offsets [mm]: "
print ["{0:0.3f}".format(i) for i in offsets]


# Quickly open the PEDe file and count lines only:
lineN= sum(1 for line in open(file))          

print "Parameters from Simulation and PEDE:"
print "PEDE Trials ",lineN  
print "Total number of variable parameters", len(expectPars)
print " "

trackN = [] # track count correspond to line number 

data = [[[0 for i_data in xrange(3)] for i_lin in xrange(lineN)] for i_par in xrange(len(expectPars))]

#Open FoM file 
with open(file) as f:
	#For each line 
	i_line = 0
	i_count = 0
	for line in f:  #Line is a string
		number_str = line.split()
		#Loop over expected parameters and store
		#Always 3 element spacing (hence hard-coded 3)
		for i_par in range(0, (len(number_str)-1)/3 ):
			label=int(number_str[0+i_par*3])
			#print "label", label
			misal=float(number_str[1+i_par*3])
			error=float(number_str[2+i_par*3])
			if (label in expectPars):
				#print "label=", label, "misal=", misal, "error=", error, "i_count=", i_count, "i_line=", i_line
				data[i_count][i_line][0] = label
				data[i_count][i_line][1] = misal
				data[i_count][i_line][2] = error
				i_count +=1 
		
		trackN.append(int(number_str[-1]))
		i_line += 1
		i_count = 0
		#print i_line

#print data
# print trackN	

##################PLOTING##############################
NewOffsets=[] # new offsets 
dMData=[] # store all dM
dMPar=[] # store corresponding par
errors=[]

misX= [[0 for i_par in xrange(len(expectPars)/int(globalN))] for i_lin in xrange(lineN)]
recoX=[[0 for i_par in xrange(len(expectPars)/int(globalN))] for i_lin in xrange(lineN)]
recoXError=[[0 for i_par in xrange(len(expectPars)/int(globalN))] for i_lin in xrange(lineN)]
misY=[[0 for i_par in xrange(len(expectPars)/int(globalN))] for i_lin in xrange(lineN)]
recoY=[[0 for i_par in xrange(len(expectPars)/int(globalN))] for i_lin in xrange(lineN)]
recoYError=[[0 for i_par in xrange(len(expectPars)/int(globalN))] for i_lin in xrange(lineN)]
dMY=[[0 for i_par in xrange(len(expectPars)/int(globalN))] for i_lin in xrange(lineN)]
dMX=[[0 for i_par in xrange(len(expectPars)/int(globalN))] for i_lin in xrange(lineN)]

plt.rcParams.update({'font.size': 14})
#Plot difference for all modules
plt.figure(1) # FoM vs Tracks

#Loop over expected parameters
iModuleX=0
iModuleY=0
iLine=-1
for i_par in range(0, len(expectPars)):
	splitLabel = [int(x) for x in str(expectPars[i_par])]  # 0 = Module, 1= Parameter

	#put a line at  0 on a plot
	axes = plt.gca()
	axes.locator_params(nbins=4, axis='y')
	line = [[0,0], [trackN[lineN-1]+500, 0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')
	#Loop over data lines 
	for i_line in range(0, lineN):
		

		dM=(data[i_par][i_line][1]-mis_C[i_par])*1e3  # mm to um rad to mrad
		dMData.append(dM)
		dMPar.append(expectPars[i_par])
		#print "i_par:", expectPars[i_par], "data[i_par][i_line][1]=", data[i_par][i_line][1], "mis_C[i_par]=", mis_C[i_par], "dM= ", (data[i_par][i_line][1]-mis_C[i_par])*1e3
		errorM=data[i_par][i_line][2]*1e3
		errors.append(errorM)
		plt.errorbar(trackN[i_line], dM, yerr=errorM, color="red") # converting 1 cm = 10'000 um
		plt.plot(trackN[i_line], dM, marker="_", color="red")
		axes.set_xlim(trackN[0]-500,trackN[lineN-1]+500)
		# axes.set_xlim(trackN[0]-500,30500)
		
		#Set label on points based on the error precision, for non-fixed modules 
		errorE = '%e' %  data[i_par][i_line][2]
		if (data[i_par][i_line][1] != 0):
			#sigDigit = int(errorE.partition('-')[2])
			label = str(round_sig(data[i_par][i_line][1])) + " mm"
			NewOffsets.append(round_sig(data[i_par][i_line][1], 4))
			axes.annotate(label, (trackN[i_line], dM), fontsize=22)
		else:
			label = "Exact/Fixed"
			NewOffsets.append(0.0)
			#sigDigit = int(errorE.partition('-')[2])
			#label = str(round_sig(data[i_par][i_line][1])) + " mm"
			#NewOffsets.append(round_sig(data[i_par][i_line][1]))
			
			axes.annotate(label, (trackN[i_line], 50), fontsize=22)
		
		
		plt.xlabel("Number of Tracks", fontsize=16)
		
		if(splitLabel[3] == 1 or splitLabel[3]==2):
			axes.set_ylim(-100, 100)
			plt.ylabel("$\Delta$ Misalignment [um]", fontsize=16)
			tick_spacing = 20
			axes.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
			
		if(splitLabel[3] == 3 or splitLabel[3]==4 or splitLabel[3] == 5):
			axes.set_ylim(-10, 10)
			plt.ylabel("$\Delta$ Misalignment [urad]", fontsize=16)
			tick_spacing = 2
			axes.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
			
		if (splitLabel[3] == 1):
			plt.title('FoM M%s X'  %(int(splitLabel[0])) , fontsize=18)
			misX[i_line][iModuleX]=(float(mis_C[i_par])*1e3)
			#print "misX[i_line][iModuleX]= ", misX[i_line][iModuleX], "iModuleX=", iModuleX, "i_line=", i_line
			recoX[i_line][iModuleX]=(float(data[i_par][i_line][1])*1e3)
			recoXError[i_line][iModuleX]=errorM
			#print "recoXError[i_line][iModuleX]=", recoXError[i_line][iModuleX]
			dMX[i_line][iModuleX]=dM
			if (i_line == lineN-1) :
				iModuleX+=1
				if (iModuleX==8):
					iModuleX=0
			
		if(splitLabel[3]==2):
			plt.title('FoM M%s Y'  %(int(splitLabel[0])) , fontsize=18)
			misY[i_line][iModuleY]=(float(mis_C[i_par])*1e3)
			recoY[i_line][iModuleY]=(float(data[i_par][i_line][1])*1e3)
			recoYError[i_line][iModuleY]=errorM
			dMY[i_line][iModuleX]=dM
			if (i_line == lineN-1) :
				iModuleY+=1
				if (iModuleY==8):
					iModuleY=0
		
		if(splitLabel[3]==3):
			plt.title('FoM M%s $\Phi$'  %(int(splitLabel[0])) , fontsize=18)

		if(splitLabel[3]==4):
			plt.title('FoM M%s $\Psi$'   %(int(splitLabel[0])) , fontsize=18)
		
		if(splitLabel[3]==5):
			plt.title('FoM M%s $\Theta$'   %(int(splitLabel[0])) , fontsize=18)

		
		

	#plt.savefig(str(expectPars[i_par])+".png")
	plt.clf()

recoXMean=[]
recoXMeanError=[]
recoYMean=[]
recoYMeanError=[]
dMXMean=[]
dMYMean=[]

offsetsX = [0 for i in xrange(0,8)]
offsetsY = [0 for i in xrange(0,8)]


if (globalN == 1):
	offsetsX=offsets
if (globalN==2):
	offsetsX=offsets[0::2]
	offsetsY=offsets[1::2]


for i_module in range(0, 8):
	meanX=0
	meanXError=[]
	meanY=0
	meanYError=[]
	for i_line in range(0, lineN):
		meanX+=recoX[i_line][i_module]
		meanY+=recoY[i_line][i_module]
		meanXError.append(recoX[i_line][i_module])
		meanYError.append(recoY[i_line][i_module])

	newMeanX = meanX/lineN+(offsetsX[i_module]*1e3)
	recoXMean.append(newMeanX)
	print "newMeanX= ", newMeanX, " meanX/lineN",  meanX/lineN, "offsetsX[i_module]=", offsetsX[i_module]*1e3
	newMeanY = meanY/lineN+(offsetsY[i_module])*1e3
	recoYMean.append(newMeanY)
	print "newMeanY= ", newMeanY, " meanY/lineN=", meanY/lineN, "offsetsY[i_module]=", offsetsY[i_module]*1e3
	recoXMeanError.append(stats.sem(meanXError))
	recoYMeanError.append(stats.sem(meanYError))
	dMXMean.append(meanX/lineN-misX[0][i_module])
	dMYMean.append(meanY/lineN-misY[0][i_module])


colours = ["green", "blue", "black", "orange", "purple"]
spacing = [2, 3.5, 4.5, 5.5, 6.5]
yMin = -300
yMax = 300
plt.subplot(211) # X 
plt.rcParams.update({'font.size': 8})
axes = plt.gca()
axes.set_xlim(0.5, 8.5)
axes.set_ylim(yMin, yMax)
plt.title("Misalignment X", fontsize=8)
plt.ylabel("Misalignment [um]", fontsize=8)
line = [[0.5,0.0], [8.5, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
for i_module in range(0, 8):
	plt.plot(i_module+1, misX[0][i_module], marker=".", color="red")
	if(useOffsets):
		plt.plot(i_module+1, offsetsX[i_module]*1e3, marker="+", color="black")

for i_module in range(0, 8):
	line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	plt.errorbar(i_module+1, recoXMean[i_module], yerr=recoXMeanError[i_module],  color=str(colours[4]), markersize=12, elinewidth=2)
	plt.plot(i_module+1, recoXMean[i_module], marker="+", color=str(colours[4]), markersize=8)
	#print "recoXMean[i_module]=", recoXMean[i_module]
	if (recoXMeanError[i_module] == 0.0):
		plt.plot(i_module+1, recoXMean[i_module], marker="*", color=str(colours[4]), markersize=12)

#Legend (top legend)
#textstr = "Truth"
#plt.text(1, 400, textstr, color="red", fontsize=10, fontweight='bold')
if (extraLabel != -1):
	plt.text(2.0, 420, extraLabel, color="purple", fontsize=8, fontweight='bold')
plt.subplots_adjust(top=0.85)


#Legend (stats X)
avgMeanMis = sum(np.array(misX[0]))/float(len(np.array(misX[0])))
SDMis = np.std(np.array(misX[0]))
#textstr = '<Truth>=%s um\nSD Truth=%s um \n'%(int(round(avgMeanMis)), int(round(SDMis)))
#plt.text(8.7, 150, textstr, fontsize=10, color="red")
avgMeandReconTruth = sum(np.array(dMXMean))/float(len(np.array(dMXMean)))
SDdReconTruth = np.std(np.array(dMXMean))
textstrReco = "<(Mean - Tr.>={0}".format(int(round(avgMeandReconTruth)))+ "um\nSD (Mean - Tr.)={0}".format(int(round(SDdReconTruth)))+" um \n" 
plt.text(8.6, 50-100, textstrReco, fontsize=8, color=str(colours[4]))
plt.subplots_adjust(right=0.77)

plt.xticks(fontsize=8, rotation=0)
plt.yticks(fontsize=8, rotation=0)
plt.xlabel("Module", fontsize=8)

plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)

plt.subplot(212) # Y 
plt.rcParams.update({'font.size': 8})
axes2 = plt.gca()
axes2.set_xlim(0.5, 8.5)
axes2.set_ylim(yMin, yMax)
plt.title("Misalignment Y", fontsize=8)
plt.ylabel("Misalignment [um]", fontsize=8)
line = [[0.5,0.0], [8.5, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
for i_module in range(0, 8):
	plt.plot(i_module+1, misY[0][i_module], marker=".", color="red")
	if(useOffsets):
		plt.plot(i_module+1, offsetsY[i_module]*1e3, marker="+", color="black")

for i_module in range(0, 8):
	line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	plt.errorbar(i_module+1, recoYMean[i_module], yerr=recoYMeanError[i_module],  color=str(colours[4]), markersize=12, elinewidth=2)
	plt.plot(i_module+1, recoYMean[i_module],  marker="+", color=str(colours[4]), markersize=8)
	if (recoYMeanError[i_module] == 0.0):
		plt.plot(i_module+1, recoYMean[i_module], marker="*", color=str(colours[4]), markersize=12)

#Legend (stats Y)
avgMeanMis = sum(np.array(misY[0]))/float(len(np.array(misY[0])))
SDMis = np.std(np.array(misY[0]))
#textstr = '<Truth>=%s um\nSD Truth=%s um \n'%(int(round(avgMeanMis)), int(round(SDMis)))
#plt.text(8.7, 150, textstr, fontsize=10, color="red")
avgMeandReconTruth = sum(np.array(dMYMean))/float(len(np.array(dMYMean)))
SDdReconTruth = np.std(np.array(dMYMean))
textstrReco = "<(Mean - Tr.>={0}".format(int(round(avgMeandReconTruth)))+ "um\nSD (Mean - Tr.)={0}".format(int(round(SDdReconTruth)))+" um \n" 
plt.text(8.6, 50-100, textstrReco, fontsize=8, color=str(colours[4]))
plt.subplots_adjust(right=0.77)
plt.subplots_adjust(bottom=0.1)
plt.subplots_adjust(hspace=0.35)
plt.xlabel("Module", fontsize=8)


# plt.show()

if (extraLabel == -1):
	plt.savefig("XY.png")
else:
	plt.savefig(str(extraLabel)+".png")


print "Mean X:", np.array(recoXMean)*1e-3, "[mm]"
print "Mean Y:", np.array(recoYMean)*1e-3, "[mm]"

sugMean  = []
for i in range(0, len(recoXMean)):
	sugMean.append(recoXMean[i])
	sugMean.append(recoYMean[i])

offsest = " "
for i in range(0, len(sugMean)):
	offsest += str(round(sugMean[i])/1e3) + " "
	

print "New suggested Offsets from PEDE [mm]: ", offsest


print "Plots saved from:" , str(file) , "on", strftime("%Y-%m-%d %H:%M:%S")