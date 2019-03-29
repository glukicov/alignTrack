####################################################################
# Sanity plots for Tracker Alignment.  
# FoM Plots for comparison of actual misalignment vs PEDE results 
# Used as a final step of FixedFoM.py script 
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 26 November 2018 by Gleb
#####################################################################
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

#custom rounding function 
def round_sig(x, sig=2):
	return round(x, sig-int(floor(log10(abs(x))))-1)

#define input arguments 
parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--mode', help='mode', default="PEDE_Mis_art.txt")  # input file mode 
parser.add_argument('-eL', '--eL', help='label', default="None") # if extra label is passed for plot tittle 
parser.add_argument('-s', '--stationN', help='station number') 
parser.add_argument('-c', '--case', help="case study", default=None) # defined below 
parser.add_argument('-i', '--iteration', help="iteration number", default=0) # define below 
args = parser.parse_args()

print("Writing new results file!")
subprocess.call(["python" , "/Users/gleb/software/alignTrack/mpIIDESY/ConcatenatePEDE.py", "-m", "w"])

#get the passed arguments 
extraLabel=str(args.eL)
stationN = str(args.stationN)
if(extraLabel == "None"):
	extraLabel="S"+stationN
case = str(args.case)
iteration = int(args.iteration)
file = str(args.mode)


#Expected parameters for stations 
if (stationN == "0"):
	expectPars = (1011, 1012, 1021, 1022, 1031, 1032, 1041, 1042, 1051, 1052, 1061, 1062, 1071, 1072, 1081, 1082)  # XY  [Station 0]
	locationStation="lower left"
if (stationN == "12"):
	expectPars = (1211, 1212, 1221, 1222, 1231, 1232, 1241, 1242, 1251, 1252, 1261, 1262, 1271, 1272, 1281, 1282)  # XY  [Station 12]
	locationStation="upper right"
if (stationN == "18"):
	expectPars = (1811, 1812, 1821, 1822, 1831, 1832, 1841, 1842, 1851, 1852, 1861, 1862, 1871, 1872, 1881, 1882)  # XY  [Station 18]
	locationStation="lower right"

#Truth Misalignment (assume initial knowledge of misalignment)
T_mis_C=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) # No Assumed Misalignment
useOffsets=False # no offsets unless set
useTruth=False 

#Define some constants
moduleN=8
moduleArray=np.arange(1, moduleN+1) #(1, 2,...,8) for plotting 
globalN=int(len(expectPars))/int(moduleN)
offsets = [0 for i in list(range(int(moduleN*globalN)))] # first set to 0 unless doing iterations
strawModuleZPosition = [ 0.0, 134.36, 268.72, 403.08, 537.42, 671.77, 806.10, 940.406] # station CS 

#If simulation, add truth misalignment to plots 
# the yMaxY etc. (axis range) are passed here for later plotting 
## ----------------------------DS0 ---------------------
if (case =="DS0"):
	useTruth = True 
	yMaxY =  120
	yMinY = -120
	yMaxX =  120
	yMinX = -120
	name = "MDC1 DS1"
	if (stationN == "12"):
		T_mis_C= (0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) #S12
	if (stationN == "18"):
		T_mis_C= (0.0, 0.0, 0.059, 0.059, 0.004, 0.007, 0.033, -0.044, -0.064, 0.066, 0.095, 0.028, -0.071, -0.023, -0.037, -0.04) #S18
	
	if(iteration == 2):
		useOffsets = True
		name="MDC1 DS0 Iteration 2"
		if (stationN == "12"):
			offsets=( -0.019, -0.007, -0.018, 0.006, 0.061, 0.001, -0.013, -0.001, -0.005, 0.006, 0.001, 0.0, 0.0, -0.003, -0.008, -0.002) #S12
		if (stationN == "18"):
			offsets=( 0.012, -0.025, 0.028, 0.026, -0.028, -0.001, -0.006, -0.04, -0.068, 0.05, 0.074, 0.022, -0.032, -0.007, 0.02, -0.024) #S18

	if(iteration == 3):
		useOffsets = True
		name="MDC1 DS0 Iteration 3"
		if (stationN == "12"):
			offsets=( -0.018, 0.009, -0.023, 0.007, 0.072, -0.01, -0.018, -0.006, -0.009, -0.008, -0.001, 0.0, 0.002, 0.0, -0.006, 0.009) #S12
		if (stationN == "18"):
			offsets=( -0.005, -0.023, 0.033, 0.019, -0.025, -0.007, 0.006, -0.038, -0.073, 0.062, 0.097, 0.032, -0.042, -0.016, 0.008, -0.029) #S18
# ##--------------------------------------------

## ----------------------------DS1---------------------------
if (case =="DS1"):
	useTruth = True 
	name="MDC1 DS1"
	if (stationN == "12"):
		T_mis_C= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) # No Assumed Misalignment
		yMaxY =  120
		yMinY = -120
		yMaxX =  120
		yMinX = -120
	if (stationN == "18"):
		T_mis_C= (0.082, 0.013, -0.102, -0.01, -0.192, 0.138, -0.192, 0.148, -0.264, 0.129, -0.168, 0.188, -0.002, 0.081, 0.03, -0.097) #S18
		yMaxY =  320
		yMinY = -320
		yMaxX =  320
		yMinX = -320
## --------------------------------------------------------

#### DATA offsets only 
## ----------------------------15922---------------------------
if (case =="15922"):
	yMaxY =  220
	yMinY = -220
	yMaxX =  120
	yMinX = -120
	name = " "
# ##----------------------------------------

if ( len(offsets) != len(expectPars) ):
	print("Enter parameter data in the right format!")

#Now print and check for correctly loaded offsets 
print("Station S:",stationN)
print("Expected Parameters: ", expectPars)
print("Using truth:", useTruth)
print("Using offsets:", useOffsets)
# input("Correct? [press enter]") 
print("Truth Misalignment ", name, stationN, ": ", T_mis_C)
# input("Truth Misalignment correct? [press enter]")
if (useOffsets==True):
	print("Offsets [mm] ", name , offsets)
	# input("Offsets correct? [press enter]")

 # the truth is the only known misalignment 
mis_C=T_mis_C 

# Quickly open the PEDE file and count lines only:
lineN= sum(1 for line in open(file))          

print("Parameters from Simulation and PEDE:")
print("PEDE Trials ",lineN  )
print("Total number of variable parameters", len(expectPars))
print(" ")

#Containers 
trackN = [] # track count correspond to line number 
data = [[[0 for i_data in range(3)] for i_lin in range(lineN)] for i_par in range(len(expectPars))]

#Open FoM file 
with open(file) as f:
	#For each line 
	i_line = 0
	i_count = 0
	for line in f:  #Line is a string
		number_str = line.split()
		#Loop over expected parameters and store
		#Always 3 element spacing (hence hard-coded 3)
		for i_par in range(0, int((len(number_str)-1)/3)  ):
			label=int(number_str[0+i_par*3])
			misal=float(number_str[1+i_par*3])
			error=float(number_str[2+i_par*3])
			if (label in expectPars):
				data[i_count][i_line][0] = label
				data[i_count][i_line][1] = misal
				data[i_count][i_line][2] = error
				i_count +=1 
		#per line add counters 
		trackN.append(int(number_str[-1]))
		i_line += 1
		i_count = 0

##################PLOTING##############################
NewOffsets=[] # new offsets 
dMData=[] # store all dM
dMPar=[] # store corresponding par
errors=[] # parameter errors 

# Containers 
misX= [[0 for i_par in range(int(len(expectPars)/int(globalN)))] for i_lin in range(lineN)]
recoX=[[0 for i_par in range(int(len(expectPars)/int(globalN)))] for i_lin in range(lineN)]
recoXError=[[0 for i_par in range(int(len(expectPars)/int(globalN)))] for i_lin in range(lineN)]
misY=[[0 for i_par in range(int(len(expectPars)/int(globalN)))] for i_lin in range(lineN)]
recoY=[[0 for i_par in range(int(len(expectPars)/int(globalN)))] for i_lin in range(lineN)]
recoYError=[[0 for i_par in range(int(len(expectPars)/int(globalN)))] for i_lin in range(lineN)]
dMY=[[0 for i_par in range(int(len(expectPars)/int(globalN)))] for i_lin in range(lineN)]
dMX=[[0 for i_par in range(int(len(expectPars)/int(globalN)))] for i_lin in range(lineN)]

#Global fontsize for plots 
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
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	#Loop over data lines 
	for i_line in range(0, lineN):
		dM=(data[i_par][i_line][1]-mis_C[i_par])*1e3  # mm to um rad to mrad
		dMData.append(dM)
		dMPar.append(expectPars[i_par])
		errorM=data[i_par][i_line][2]*1e3
		errors.append(errorM)
		plt.errorbar(trackN[i_line], dM, yerr=errorM, color="red") # converting 1 cm = 10'000 um
		plt.plot(trackN[i_line], dM, marker="_", color="red")
		axes.set_xlim(trackN[0]-500,trackN[lineN-1]+500)
		
		#Set label on points based on the error precision, for non-fixed modules 
		errorE = '%e' %  data[i_par][i_line][2]
		if (data[i_par][i_line][1] != 0):
			label = str(round_sig(data[i_par][i_line][1])) + " mm"
			NewOffsets.append(round_sig(data[i_par][i_line][1], 4))
			axes.annotate(label, (trackN[i_line], dM), fontsize=22)
		else:
			label = "Exact/Fixed"
			NewOffsets.append(0.0)
			axes.annotate(label, (trackN[i_line], 50), fontsize=22)
		
		plt.xlabel("Number of Tracks", fontsize=16)
		
		#Now set y-label according to parameter type 
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
			
		#Now set title according to parameter type 
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

	#save of close the figure 
	#plt.savefig(str(expectPars[i_par])+".png")
	plt.clf()

#New mean containers 
recoXMean=[]
recoXMeanError=[]
recoYMean=[]
recoYMeanError=[]
dMXMean=[]
dMYMean=[]
offsetsX = [0 for i in list(range(0,8))]
offsetsY = [0 for i in list(range(0,8))]

#Check number of parameters 
if (globalN == 1):
	offsetsX=offsets
if (globalN==2):
	offsetsX=offsets[0::2]
	offsetsY=offsets[1::2]
if (globalN>2):
	input("Not ready for so many parameters - edit the script! [press enter]") 

#Calculate the new suggested parameter, based on the mean 
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
	print("newMeanX= ", newMeanX, " meanX/lineN",  meanX/lineN, "offsetsX[i_module]=", offsetsX[i_module]*1e3)
	newMeanY = meanY/lineN+(offsetsY[i_module])*1e3
	recoYMean.append(newMeanY)
	print("newMeanY= ", newMeanY, " meanY/lineN=", meanY/lineN, "offsetsY[i_module]=", offsetsY[i_module]*1e3)
	recoXMeanError.append(stats.sem(meanXError))
	recoYMeanError.append(stats.sem(meanYError))
	dMXMean.append(meanX/lineN-misX[0][i_module])
	dMYMean.append(meanY/lineN-misY[0][i_module])


#################################
### Correct the angle in truth misalignment

# TODO implement in the loop + plotting 
#form X (z-position), Y (misalignment) array 
dataRad=[strawModuleZPosition, misX[0]] 
dataVer=[strawModuleZPosition, misY[0]]
dataArray = (dataRad, dataVer)
MisCorrected=[]

print("********Uncorrected Truth******")
sys.stdout.write("services.Geometry.strawtracker.strawModuleRShift"+str(stationN)+": [ ")
for item in misX[0][:-1]:
	sys.stdout.write( str(round(int(item)*1e-3,3)) + ", ")  # um -> mm 
sys.stdout.write( str( round(int(misX[0][-1])*1e-3,3) ) )
sys.stdout.write(" ]\n")

sys.stdout.write("services.Geometry.strawtracker.strawModuleHShift"+str(stationN)+": [ ")
for item in misY[0][:-1]:
    sys.stdout.write( str( round(int(item)*1e-3,3)) + ", " ) 
sys.stdout.write( str( round(int(misY[0][-1])*1e-3,3) ) ) 
sys.stdout.write(" ]\n")
sys.stdout.flush()
print("********************")

for i_global in range(0, int(globalN)):
	
	data = dataArray[i_global]
	#Do a linear fit 
	x_new = np.linspace(float(min(data[0])), float(max(data[0])), num=1e3) # generate x-points for evaluation 
	coefs = poly.polyfit(data[0], data[1], 1) # straight line
	ffit = poly.polyval(x_new, coefs) # plot over generated points 

	#Correct the truth
	corrected_T_mis_C = []
	for i_module in range(0, len(moduleArray)):
	    # corr = Y + Gradient * X + Intercept 
	    correction = data[1][i_module] - coefs[1] * data[0][i_module] - coefs[0]
	    corrected_T_mis_C.append(correction)
	
	MisCorrected.append(corrected_T_mis_C)



##################################

colours = ["green", "blue", "black", "orange", "purple"]
plt.subplot(211) # X 
plt.rcParams.update({'font.size': 8})
axes = plt.gca()
axes.set_xlim(0.5, 8.5)
axes.set_ylim(yMinX, yMaxX)
xTitle = extraLabel + " Misalignment X"
plt.title(str(xTitle), fontsize=8)
plt.ylabel("Misalignment [um]", fontsize=8)
plt.minorticks_on()
axes.tick_params(axis='x',which='minor',bottom=False)
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
line = [[0.5,0.0], [8.5, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')

plt.plot(moduleArray, np.array(MisCorrected[0]), marker=".", color="red", label="Truth Mis.")
if(useOffsets):
	plt.plot(moduleArray, np.array(offsetsX)*1e3, marker="+", color="black", label="Reco. Mis.\n(previous iteration)")

for i_module in range(0, 8):
	line = [[i_module+0.5,yMinX], [i_module+0.5, yMaxX]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	print("recoXMean[i_module]", recoXMean[i_module], "recoXMeanError[i_module]", recoXMeanError[i_module])
	if (recoXMean[i_module] == 0.0):
		plt.plot(i_module+1, 0.0, marker="*", color=str(colours[4]), markersize=12) # if module was fixed - plot a star 
print(moduleArray)
print(recoXMean)
plt.errorbar(moduleArray, np.array(recoXMean), yerr=recoXMeanError,  color=str(colours[4]), markersize=12, elinewidth=2, label="Reco. Mis.\n(this iteration)")
plt.plot(moduleArray, np.array(recoXMean), marker="+", color=str(colours[4]))

#Legend (stats X)
if(useTruth == True):
	avgMeanMisX = sum(np.abs(np.array(MisCorrected[0])))/float(len(np.array(MisCorrected[0])))
	textstr = '<Truth X>=%s um'%(int(round(avgMeanMisX)))
	plt.text(8.7, yMaxX, textstr, fontsize=10, color="red")
	avgMeanRecoX = sum(np.abs(np.array(recoXMean)))/float(len(np.array(recoXMean)))
	textstr = '<Reco X>=%s um'%(int(round(avgMeanRecoX)))
	plt.text(8.7, yMaxX*0.8, textstr, fontsize=10, color="purple")

# plt.legend(loc=locationStation)
axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xticks(fontsize=8, rotation=0)
plt.yticks(fontsize=8, rotation=0)
plt.xlabel("Module", fontsize=8)

plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)

plt.subplot(212) # Y 
plt.rcParams.update({'font.size': 8})
axes2 = plt.gca()
axes2.set_xlim(0.5, 8.5)
axes2.set_ylim(yMinY, yMaxY)
yTitle = extraLabel + " Misalignment Y"
plt.title(str(yTitle), fontsize=8)
plt.ylabel("Misalignment [um]", fontsize=8)
plt.minorticks_on()
axes2.tick_params(axis='x',which='minor',bottom=False)
axes2.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
line = [[0.5,0.0], [8.5, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')

plt.plot(moduleArray, np.array(MisCorrected[1]), marker=".", color="red", label="Truth Mis.")
if(useOffsets):
	plt.plot(moduleArray, np.array(offsetsY)*1e3, marker="+", color="black", label="Reco. Mis.\n(previous iteration)")

for i_module in range(0, 8):
	line = [[i_module+0.5,yMinY], [i_module+0.5, yMaxY]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	print("recoYMean[i_module]", recoYMean[i_module], "recoYMeanError[i_module]", recoYMeanError[i_module])
	if (recoYMean[i_module] == 0.0):
		plt.plot(i_module+1, 0.0, marker="*", color=str(colours[4]), markersize=12)

plt.errorbar(moduleArray, np.array(recoYMean), yerr=recoYMeanError[i_module],  color=str(colours[4]), markersize=12, elinewidth=2, label="Reco. Mis.\n(this iteration)")
plt.plot(moduleArray, np.array(recoYMean),  marker="+", color=str(colours[4]), markersize=8)

#Legend (stats Y)
if(useTruth):
	avgMeanMisY = sum(np.abs(np.array(MisCorrected[1])))/float(len(np.array(misY[0])))
	textstr = '<Truth Y>=%s um'%(int(round(avgMeanMisY)))
	plt.text(8.7, yMaxY, textstr, fontsize=10, color="red")
	avgMeanRecoY = sum(np.abs(np.array(recoYMean)))/float(len(np.array(recoYMean)))
	textstr = '<Reco Y>=%s um'%(int(round(avgMeanRecoY)))
	plt.text(8.7, yMaxY*0.8, textstr, fontsize=10, color="purple")
	
axes2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.subplots_adjust(bottom=0.1)
plt.subplots_adjust(hspace=0.35)
plt.xlabel("Module", fontsize=8)


if (extraLabel == -1):
	plt.savefig("XY.png")
else:
	plt.savefig(str(extraLabel)+".png", dpi=200)

sugMean  = []
#Take away the 0th module shift from the rest 
for i in range(0, len(recoXMean)):
	sugMean.append(recoXMean[i])
	sugMean.append(recoYMean[i])

offset = " "
for i in range(0, (len(sugMean)-1) ):
	offset += str(round(sugMean[i])/1e3) + ", "
offset += str(round(sugMean[-1])/1e3)
	

print("New suggested Offsets from PEDE [mm]: ", offset)

#Arrange offsets from X1, Y1, X2....X8, Y8 into strawModuleRShift12, etc. 

recoXMean = np.array(recoXMean)*1e-3
recoYMean = np.array(recoYMean)*1e-3
 
radialOffsetPerModule=np.round_(recoXMean, 3) #X per station
verticalOffsetPerModule=np.round_(recoYMean, 3) #Y per station

print("radialOffsetPerModule", radialOffsetPerModule)
print("verticalOffsetPerModule", verticalOffsetPerModule)

f=open("OffsetsPerModule"+str(stationN)+".fcl", "w+")

f.write("//StrawOffsets S"+str(stationN)+"\n")
f.write("//(X, Y):\n")
f.write("// offsets=("+str(offset)+") #S"+str(stationN)+"\n")


f.write("services.Geometry.strawtracker.strawModuleRShift"+str(stationN)+": [ ")
for item in radialOffsetPerModule[:-1]:
	f.write( str(item) + ", ")  # um -> mm 
f.write( str( radialOffsetPerModule[-1] ) )
f.write(" ]\n")

f.write("services.Geometry.strawtracker.strawModuleHShift"+str(stationN)+": [ ")
for item in verticalOffsetPerModule[:-1]:
    f.write( str( item) + ", " ) 
f.write( str( verticalOffsetPerModule[-1] ) ) 
f.write(" ]\n")

f.close()

print("Offsets written to file OffsetsPerModule"+str(stationN)+".fcl")
print("Plots saved from:" , str(file) , "on", strftime("%Y-%m-%d %H:%M:%S"))