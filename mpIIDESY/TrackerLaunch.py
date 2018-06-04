#!/usr/bin/python

####################################################################
# Sanity plots for Tracker Alignment.  
# FoM Plots for comparison of actual misalignment vs PEDE results 
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 8 May 2018 by Gleb
#####################################################################
import argparse, sys

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--mode', help='mode')
args = parser.parse_args()
from math import log10, floor

def round_sig(x, sig=2):
	return round(x, sig-int(floor(log10(abs(x))))-1)

file = str(args.mode)

import matplotlib.pyplot as plt #for plotting 
import matplotlib.ticker as ticker
import numpy as np  # smart arrays 
import itertools # smart lines
from time import gmtime, strftime 
import subprocess


# #Truth Misalignment 
# expectPars = (11, 12, 21, 22, 31, 32, 51, 52, 71, 72, 81, 82)
# mis_C = (0.1, 0.15, 0.05, 0.05, -0.1, -0.15, -0.07, 0.1, 0.05, 0.07, 0.0, 0.0)

expectPars = (11, 12, 21, 22, 31, 32, 41, 42, 51, 52, 61, 62, 71, 72, 81, 82)

mis_C=(-0.2, 0.1, 0.08, 0.15, 0.2, -0.1, -0.25, 0.3, 0.15, 0.2, 0.1, -0.25, 0.2, 0.07, -0.06, 0.06)

# T_mis_C = (0.1, 0.15, 0.05, 0.05, -0.1, -0.15, 0.0, 0.0, -0.07, 0.1, 0.0, 0.0, 0.05, 0.07, 0.0, 0.0)
# mis_C=[]
# offsets = (0.014, 0.11, 0.0, 0.0, -0.12, -0.17, 0.0, 0.0, -0.057, 0.088, 0.017, -0.0085, 0.063, 0.07, 0.0, 0.0)
# print len(T_mis_C), len(offsets)
# for i in range(0, len(T_mis_C)):
# 	mis_C.append(float(T_mis_C[i] - offsets[i]))

print "Truth Misalignments after offsets", mis_C

# expectPars = (11, 12, 21, 22, 31, 32, 51, 52, 71, 72, 81, 82)
# mis_C = (0.1, 0.15, 0.05, 0.05, -0.1, -0.15, -0.07, 0.1, 0.05, 0.07, 0.0, 0.0)

# expectPars = (11, 12, 21, 22, 31, 32, 51, 52, 61, 62 , 71, 72, 81, 82 )
# mis_C = (0.1, 0.15, 0.05, 0.05, -0.1, -0.15, -0.07, 0.1, 0.0, 0.0, 0.05, 0.07, 0.0, 0.0)

# mis_C = (0.2, 0.15, -0.1, -0.1)
# expectPars = (21, 22, 31, 32)

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
offsests=[]

plt.rcParams.update({'font.size': 14})
#Plot difference for all modules
plt.figure(1)

#Loop over expected parameters
for i_par in range(0, len(expectPars)):

	splitLabel = [int(x) for x in str(expectPars[i_par])]  # 0 = Module, 1= Parameter

	#put a line at  0 on a plot
	axes = plt.gca()
	axes.locator_params(nbins=4, axis='y')
	line = [[0,0], [trackN[lineN-1]+500, 0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green', linestyle=":")
	#Loop over data lines 

	for i_line in range(0, lineN):
		
		dM=(data[i_par][i_line][1]-mis_C[i_par])*1e3  # mm to um rad to mrad 
		#print "data[i_par][i_line][1]=", data[i_par][i_line][1], "mis_C[i_par]=", mis_C[i_par], "dM= ", (data[i_par][i_line][1]-mis_C[i_par])*1e3
		errorM=data[i_par][i_line][2]*1e3
		plt.errorbar(trackN[i_line], dM, yerr=errorM, color="red") # converting 1 cm = 10'000 um
		plt.plot(trackN[i_line], dM, marker="_", color="red")
		axes.set_xlim(trackN[0]-500,trackN[lineN-1]+500)
		# axes.set_xlim(trackN[0]-500,30500)
		
		#Set label on points based on the error precision, for non-fixed modules 
		errorE = '%e' %  data[i_par][i_line][2]
		if (data[i_par][i_line][1] != 0):
			sigDigit = int(errorE.partition('-')[2])
			label = str(round_sig(data[i_par][i_line][1])) + " mm"
			offsests.append(round_sig(data[i_par][i_line][1]))
			axes.annotate(label, (trackN[i_line], dM), fontsize=22)
		else:
			label = "Exact/Fixed"
			offsests.append(0.0)
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
			
		if(splitLabel[1]==2):
			plt.title('FoM M%s Y'  %(int(splitLabel[0])) , fontsize=18)
		
		if(splitLabel[1]==3):
			plt.title('FoM M%s $\Phi$'  %(int(splitLabel[0])) , fontsize=18)

		if(splitLabel[1]==4):
			plt.title('FoM M%s $\Psi$'   %(int(splitLabel[0])) , fontsize=18)
		
		if(splitLabel[1]==5):
			plt.title('FoM M%s $\Theta$'   %(int(splitLabel[0])) , fontsize=18)

		
		

	plt.savefig(str(expectPars[i_par])+".png")
	plt.clf()

#Now combine produced plots into a single file:
#convert -append 11.png 12.png 1.png
even=expectPars[1::2]
odd=expectPars[0::2]

newOdd=[]
newEven=[]

for i_par in range(0, int(len(expectPars)/2)):

	f1=str(odd[i_par])+".png" 
	f2=str(even[i_par])+".png"
	f3=str(i_par)+".png"
	subprocess.call(["convert" , "-append", str(f1), str(f2), str(f3)])
	if (i_par % 2 == 0):
		newEven.append(f3)
	if (i_par %2 != 0):
		newOdd.append(f3)

subprocess.call(["convert" , "+append", str(newEven[0]), str(newOdd[0]), str(newEven[1]), str(newOdd[1]),  "Row1.png"])
subprocess.call(["convert" , "+append", str(newEven[2]), str(newOdd[2]), str(newEven[3]), str(newOdd[3]),   "Row2.png"])
subprocess.call(["convert" , "-append", "Row1.png", "Row2.png", "FoM.png"])


print "Plots saved from:" , str(file) , "on", strftime("%Y-%m-%d %H:%M:%S")

offsest = " "
for i in range(0, len(offsests)):
	offsest += str(offsests[i]) + " "

print "Suggested Offsets from PEDE: ", offsest