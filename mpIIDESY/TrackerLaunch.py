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

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--mode', help='mode')
parser.add_argument('-eL', '--eL', help='mode')
args = parser.parse_args()
from math import log10, floor, ceil

def round_sig(x, sig=2):
	return round(x, sig-int(floor(log10(abs(x))))-1)

def float_round(num, places = 0, direction = floor):
    return direction(num * (10**places)) / float(10**places)

extraLabel =-1
extraLabel=str(args.eL)

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


# T_mis_C=(0.1, -0.07, -0.08, 0.05, 0.15, 0.1, -0.04, 0.01) # Case C (Initial)
# expectPars = (11, 21, 31, 41, 51, 61, 71, 81)

expectPars = (11, 12, 21, 22, 31, 32, 41, 42, 51, 52, 61, 62, 71, 72, 81, 82)

#Truth Misalignment 

T_mis_C = (0.1, 0.15, 0.05, 0.05, -0.1, -0.15, 0.0, 0.0, -0.07, 0.1, 0.0, 0.0, 0.05, 0.07, 0.0, 0.0) # Case A (Initial)

# T_mis_C=(-0.2, 0.1, 0.08, 0.15, 0.2, -0.1, -0.25, 0.3, 0.15, 0.2, 0.1, -0.25, 0.2, 0.07, -0.06, 0.06) # Case B (Initial)

# T_mis_C=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) # No Misalignment

#T_mis_C=(0.103, 0.195, 0.080, 0.150, 0.060, 0.148, 0.035, 0.121, 0.016,0.106, -0.006, 0.082, -0.028, 0.068, -0.060, 0.060 ) # Case D : Residual Mis. 

globalN=int(len(expectPars))/int(8)
offsets = [0 for i in xrange(8*globalN)]

# Run 0 
mis_C=T_mis_C  # the truth is the only misalignment 
print "Initial Truth Misalignment [mm]: ", mis_C
print "With expected Parameters: ", expectPars
# raw_input("Truth Misalignment correct? [press enter]") 

# # ----------------------------
#Run 1
#mis_C=[] # set temp to 0 

# offsets = (0.028, 0.049, -0.033, -0.041, -0.184, -0.224, -0.08, -0.075, -0.139, 0.043, -0.053, -0.05, 0.016, 0.026, -0.008, -0.031) # Mean M82F Case A (Run 1)

# offsets = (-0.146, -0.029, 0.112, 0.025, 0.207, -0.212, -0.261, 0.189, 0.116, 0.098, 0.047, -0.335, 0.129, -0.002, -0.138, 0.003) # Mean M8A Case B (Run 1)

# offsets = (0.06, -0.062, -0.035, 0.117, 0.22, 0.159, -0.002, 0.015) # Mean XY  Case C (Run 1)

#print "Offsets Run 1 [mm]: ", offsets
#raw_input("Offsets :: Run 1 correct? [press enter]") 
# for i in range(0, len(T_mis_C)):
# 	mis_C.append(float(T_mis_C[i] - offsets[i]))
#----------------------------

# ##----------------------------
# # Run 2
# T_mis_C=mis_C #set truth as the previous misalignment
# mis_C = [] 

offsets =(0.018, 0.052, -0.034, -0.034, -0.183, -0.23, -0.077, -0.09, -0.136, 0.03, -0.05, -0.064, 0.021, 0.024, -0.003, -0.029) # M82F Mean Run 2

# offsets = (-0.165, -0.031, 0.118, 0.026, 0.235, -0.208, -0.226, 0.198, 0.154, 0.105, 0.08, -0.331, 0.148, -0.0, -0.148, 0.001)  # M8A Mean Run 2


# print "Offsets Run 2 [mm]: ", offsets
# raw_input("Offsets :: Run 2 correct? [press enter]") 
# for i in range(0, len(T_mis_C)):
# 	mis_C.append(float(T_mis_C[i] - offsets[i]))
# ##----------------------------


# #----------------------------
# # Run 3
# T_mis_C=mis_C #set truth as the previous misalignment
# mis_C = [] 



# print "Offsets Run 3 [mm]: ", offsets
# raw_input("Offsets :: Run 2 correct? [press enter]") 
# for i in range(0, len(T_mis_C)):
# 	mis_C.append(float(T_mis_C[i] - offsets[i]))
# #----------------------------

print "Truth Misalignments after offsets [mm]: "
print ["{0:0.3f}".format(i) for i in mis_C]

if ( len(mis_C) != len(expectPars) ):
	print "Enter Truth data in the right format!"


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
offsests=[] # new offsets 
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
			offsests.append(round_sig(data[i_par][i_line][1], 4))
			axes.annotate(label, (trackN[i_line], dM), fontsize=22)
		else:
			label = "Exact/Fixed"
			offsests.append(0.0)
			#sigDigit = int(errorE.partition('-')[2])
			#label = str(round_sig(data[i_par][i_line][1])) + " mm"
			#offsests.append(round_sig(data[i_par][i_line][1]))
			
			axes.annotate(label, (trackN[i_line], 50), fontsize=22)
		
		
		plt.xlabel("Number of Tracks", fontsize=16)
		
		if(splitLabel[1] == 1 or splitLabel[1]==2):
			axes.set_ylim(-100, 100)
			plt.ylabel("$\Delta$ Misalignment [um]", fontsize=16)
			tick_spacing = 20
			axes.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
			
		if(splitLabel[1] == 3 or splitLabel[1]==4 or splitLabel[1] == 5):
			axes.set_ylim(-10, 10)
			plt.ylabel("$\Delta$ Misalignment [urad]", fontsize=16)
			tick_spacing = 2
			axes.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
			
		if (splitLabel[1] == 1):
			plt.title('FoM M%s X'  %(int(splitLabel[0])) , fontsize=18)
			misX[i_line][iModuleX]=(float(mis_C[i_par])*1e3)
			#print "misX[i_line][iModuleX]= ", misX[i_line][iModuleX], "iModuleX=", iModuleX, "i_line=", i_line
			if (i_line != 0):
				recoX[i_line][iModuleX]=(float(data[i_par][i_line][1])*1e3+recoX[i_line-1][iModuleX])
			else:
				recoX[i_line][iModuleX]=(float(data[i_par][i_line][1])*1e3)
			recoXError[i_line][iModuleX]=errorM
			#print "recoXError[i_line][iModuleX]=", recoXError[i_line][iModuleX]
			dMX[i_line][iModuleX]=dM
			if (i_line == lineN-1) :
				iModuleX+=1
				if (iModuleX==8):
					iModuleX=0
			
		if(splitLabel[1]==2):
			plt.title('FoM M%s Y'  %(int(splitLabel[0])) , fontsize=18)
			misY[i_line][iModuleY]=(float(mis_C[i_par])*1e3)
			if (i_line != 0):
				recoY[i_line][iModuleY]=(float(data[i_par][i_line][1])*1e3+recoY[i_line-1][iModuleY])
			else:
				recoY[i_line][iModuleY]=(float(data[i_par][i_line][1])*1e3)
			recoYError[i_line][iModuleY]=errorM
			dMY[i_line][iModuleX]=dM
			if (i_line == lineN-1) :
				iModuleY+=1
				if (iModuleY==8):
					iModuleY=0
		
		if(splitLabel[1]==3):
			plt.title('FoM M%s $\Phi$'  %(int(splitLabel[0])) , fontsize=18)

		if(splitLabel[1]==4):
			plt.title('FoM M%s $\Psi$'   %(int(splitLabel[0])) , fontsize=18)
		
		if(splitLabel[1]==5):
			plt.title('FoM M%s $\Theta$'   %(int(splitLabel[0])) , fontsize=18)

		
		

	#plt.savefig(str(expectPars[i_par])+".png")
	plt.clf()

#print recoXError
offsetsX = [0 for i in xrange(0,8)]
offsetsY = [0 for i in xrange(0,8)]

if (globalN == 1):
	offsetsX=offsets
if (globalN==2):
	offsetsX=offsets[0::2]
	offsetsY=offsets[1::2]



colours = ["green", "blue", "black", "orange", "purple"]
spacing = [2, 3.5, 4.5, 5.5, 6.5]
yMin = -300
yMax = 300
plt.subplot(211) # X 
plt.rcParams.update({'font.size': 10})
axes = plt.gca()
axes.set_xlim(0.5, 8.5)
axes.set_ylim(yMin, yMax)
plt.title("Misalignment X", fontsize=10)
plt.ylabel("Misalignment [um]", fontsize=10)
line = [[0.5,0.0], [8.5, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
for i_module in range(0, 8):
	plt.plot(i_module+1, misX[0][i_module], marker=".", color="red")
	plt.plot(i_module+1, offsetsX[i_module]*1e3, marker="+", color="black")

for i_line in range(0, lineN):
	for i_module in range(0, 8):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
		#print "recoXError[i_line][i_module]=", recoXError[i_line][i_module]
		plt.errorbar(i_module+1, recoX[i_line][i_module]+offsetsX[i_module]*1e3, yerr=recoXError[i_line][i_module],  color="purple", markersize=12, elinewidth=2)
		if (recoXError[i_line][i_module] == 0.0):
			plt.plot(i_module+1, recoX[i_line][i_module], marker="*", color=str(colours[i_line]), markersize=12)

#Legend (top legend)
textstr = "Truth"
plt.text(1, 400, textstr, color="red", fontsize=10, fontweight='bold')
if (extraLabel != -1):
	plt.text(6, 350, extraLabel, color="purple", fontsize=12, fontweight='bold')
for i_line in range(0, lineN):
	textstr = "Iteration " + str(i_line)
	plt.text(spacing[i_line], 400, textstr, color=str(colours[i_line]), fontsize=10, fontweight='bold')
plt.subplots_adjust(top=0.85)

#Legend (stats X)
avgMeanMis = sum(np.array(misX[0]))/float(len(np.array(misX[0])))
SDMis = np.std(np.array(misX[0]))
textstr = '<Truth>=%s um\nSD Truth=%s um \n'%(int(round(avgMeanMis)), int(round(SDMis)))
plt.text(8.7, 150, textstr, fontsize=10, color="red")
for i_line in range(0, lineN):
	avgMeandReconTruth = sum(np.array(dMX[i_line]))/float(len(np.array(dMX[i_line])))
	SDdReconTruth = np.std(np.array(dMX[i_line]))
	textstrReco = "<(It. "+str(i_line)+" - Tr.>={0}".format(int(round(avgMeandReconTruth)))+ "um\nSD (It. "+str(i_line)+" - Tr.)={0}".format(int(round(SDdReconTruth)))+" um \n" 
	plt.text(8.6, 50-100*i_line, textstrReco, fontsize=10, color=str(colours[i_line]))
plt.subplots_adjust(right=0.85)

plt.subplot(212) # Y 
axes = plt.gca()
axes.set_xlim(0.5, 8.5)
axes.set_ylim(yMin, yMax)
plt.title("Misalignment Y", fontsize=10)
plt.ylabel("Misalignment [um]", fontsize=10)
line = [[0.5,0.0], [8.5, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
for i_module in range(0, 8):
	plt.plot(i_module+1, misY[i_line][i_module], marker=".", color="red")
	plt.plot(i_module+1, offsetsY[i_module]*1e3, marker="+", color="black")

for i_line in range(0, lineN):
	for i_module in range(0, 8):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
		plt.errorbar(i_module+1, recoY[i_line][i_module]+offsetsY[i_module]*1e3, yerr=recoYError[i_line][i_module],  color="purple", markersize=12, elinewidth=2)
		if (recoYError[i_line][i_module] == 0.0):
			plt.plot(i_module+1, recoY[i_line][i_module], marker="*", color=str(colours[i_line]), markersize=12)

#Legend (stats Y)
avgMeanMis = sum(np.array(misY[0]))/float(len(np.array(misY[0])))
SDMis = np.std(np.array(misY[0]))
textstr = '<Truth>=%s um\nSD Truth=%s um \n'%(int(round(avgMeanMis)), int(round(SDMis)))
plt.text(8.7, 150, textstr, fontsize=10, color="red")
for i_line in range(0, lineN):
	avgMeandReconTruth = sum(np.array(dMY[i_line]))/float(len(np.array(dMY[i_line])))
	SDdReconTruth = np.std(np.array(dMY[i_line]))
	textstrReco = "<(It. "+str(i_line)+" - Tr.>={0}".format(int(round(avgMeandReconTruth)))+ "um\nSD (It. "+str(i_line)+" - Tr.)={0}".format(int(round(SDdReconTruth)))+" um \n" 
	plt.text(8.6, 50-100*i_line, textstrReco, fontsize=10, color=str(colours[i_line]))
plt.subplots_adjust(right=0.78)

plt.xlabel("Module", fontsize=12)
if (extraLabel == -1):
	plt.savefig("XY.png")
else:
	plt.savefig(str(extraLabel)+".png")


#Now combine produced plots into a single file:
#convert -append 11.png 12.png 1.png
even=expectPars[1::2]
odd=expectPars[0::2]

newOdd=[]
newEven=[]

# for i_par in range(0, int(len(expectPars)/2)):

# 	f1=str(odd[i_par])+".png" 
# 	f2=str(even[i_par])+".png"
# 	f3=str(i_par)+".png"
# 	subprocess.call(["convert" , "-append", str(f1), str(f2), str(f3)])
# 	if (i_par % 2 == 0):
# 		newEven.append(f3)
# 	if (i_par %2 != 0):
# 		newOdd.append(f3)

# subprocess.call(["convert" , "+append", str(newEven[0]), str(newOdd[0]), str(newEven[1]), str(newOdd[1]),  "Row1.png"])
# subprocess.call(["convert" , "+append", str(newEven[2]), str(newOdd[2]), str(newEven[3]), str(newOdd[3]),   "Row2.png"])
# subprocess.call(["convert" , "-append", "Row1.png", "Row2.png", "FoM.png"])


#------DM Histo---------
#


cDM = TCanvas("cDM", "cDM", 700, 700)
cDM.Divide(2,1)
minF=min(dMData)-0.1
maxF=max(dMData)+0.1
hDMx = TH1F("hDMx", "PEDE - Truth Alignment in X; X Misalignment [um]", 49, minF, maxF )
hDMy = TH1F("hDMy", "PEDE - Truth Alignment in Y; Y Misalignment [um]", 49, minF, maxF )

# Fill histos for align-able modules (error =0 == fixed module)
for i_line in range(0, lineN):
	for i in range(0, len(dMData)):
		splitLabel = [int(x) for x in str(dMPar[i])]  # 0 = Module, 1= Parameter
		if (splitLabel[1] == 1 and errors[i] !=0):
			#print "dMPar[i]=", dMPar[i], "dMData[i]=", dMData[i]
			hDMx.Fill(dMData[i])
		if (splitLabel[1] == 2 and errors[i] !=0):
			#print "dMPar[i]=", dMPar[i], "dMData[i]=", dMData[i]
			hDMy.Fill(dMData[i])
		  	

#Set function to fit
# func = TF1("fu#c", "gaus(0)", minF, maxF)
#func->SetParameters(0.0, hMean);

#Fit function
cDM.cd(1)
hDMx.Draw() #Set errors on all bins
cDM.cd(2)
hDMy.Draw() #Set errors on all bins
#hDM.Draw("E1") #Set errors on all bins
#hDM.Fit("func")
gStyle.SetOptStat("ourRmMe"); #over/under -flows, Rms and Means with errors, number of entries
gStyle.SetOptFit(1111);  #probability, Chi2, errors, name/values of parameters
gStyle.SetStatFormat("11.4f");  # 4 sig.fig, f=float

#Save canvas as .png file
cDM.Modified()
cDM.Update()
cDM.Print("dM.png")
cDM.Print("dM.root")

xSD=hDMx.GetRMS()
xSDError=hDMx.GetRMSError()
xMean=hDMx.GetMean()
xMeanError=hDMx.GetMeanError()
ySD=hDMy.GetRMS()
ySDError=hDMy.GetRMSError()
yMean=hDMy.GetMean()
yMeanError=hDMy.GetMeanError()

print "xSD=",round(xSD),"+/-",round(xSDError), "xMean=",round(xMean),"+/-",round(xMeanError)
print "ySD=",round(ySD),"+/-",round(ySDError), "yMean=",round(yMean),"+/-",round(yMeanError)


offsest = " "
for i in range(0, len(offsests)):
	if (errors[i] != 0.0):
		offsest += str(round(offsests[i]*1e3)/1e3) + " "
	else:
		offsest += str(offsests[i]) + " "


print "New suggested Offsets from PEDE [mm]: ", offsest


# offsests_um=[]
# mis_C_um=[]
# for i in range(0, len(offsests)):
# 	#print offsests[i]
# 	#print int(round(offsests[i]*1e3))
# 	offsests_um.append(round(offsests[i]*1e3))
# 	mis_C_um.append(round(T_mis_C[i]*1e3))

# trace = go.Table(
#     header=dict(values=['Labels', 'Truth [um]', 'Iteration 0 [um]'],
#                 line = dict(color='#7D7F80'),
#                 fill = dict(color='#a1c3d1'),
#                 align = ['left'] * 5),
#     cells=dict(values=[expectPars, mis_C_um, offsests_um],
#                line = dict(color='#7D7F80'),
#                fill = dict(color='#EDFAFF'),
#                align = ['left'] * 5))

# layout = dict(width=500, height=600)
# data = [trace]
# fig = dict(data=data, layout=layout)
# plotly.plotly.iplot(fig, filename = 'styled_table')


print "Plots saved from:" , str(file) , "on", strftime("%Y-%m-%d %H:%M:%S")