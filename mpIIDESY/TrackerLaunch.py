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

file = str(args.mode)

import matplotlib.pyplot as plt #for plotting 
import matplotlib.ticker as ticker
import numpy as np  # smart arrays 
import itertools # smart lines
from time import gmtime, strftime 


#Truth Misalignment 

#Put only the active alignment parameters 
# mis_C = ( 0.2, -0.1, -0.0122, 0.0087,  -0.1, 0.15, 0.00873, -0.0070) 
# mis_C = (0.15, 0.2, -0.2, -0.1)
# mis_C = (0.08, 0.045, 0.0087, 0.0087, -0.07, -0.06, -0.0175, -0.0070) 
# mis_C = (0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0) 
# mis_C = ( -0.1, 0.15, 0.00873, -0.0070,  0.2, -0.1, -0.0122, 0.0087) 
# mis_C = ( -0.1, 0.15, 0.00873, -0.0070,  0.2, -0.1, -0.0122, 0.0087) 
# mis_C = (0.08, 0.045, 0.0087, 0.0087, -0.07, -0.06, -0.0175, -0.0070)
# mis_C = (0.1, -0.08, -0.1, 0.06) 

# mis_C = ( 0.2, -0.1, 0.0087, 0.0087,  -0.1, 0.15, -0.0175, -0.0070) 

# mis_C = (0.05, -0.1, -0.15)
# mis_C = (0.05, -0.1, -0.15)

# mis_C = ( 0.2, 0.045, -0.0122, 0.0087,  -0.1, -0.06, 0.00873, -0.0070, 0.0, 0.0, 0.0, 0.0) 

# expectPars=(11, 12, 14, 15, 41, 42, 44, 45)
# expectPars=(21, 22, 23, 25, 31, 32, 33, 35)
# expectPars=(21, 22, 24, 25, 31, 32, 34, 35)
# expectPars = (21, 22, 23, 25, 31, 32, 33, 35)
# expectPars=(31, 32, 33, 35, 41, 42, 43, 45)




mis_C = (0.1, -0.09, -0.15, 0.1)

expectPars=(21, 22, 31, 32)


constN = 1 #Number of fixed modules 
moduleN = 4 #Total number of modules

alignN = moduleN - constN
#TODO load from Sim
parN=int(len(expectPars)/alignN)# Alignment parameters  

# Quickly open the PEDe file and count lines only:
lineN= sum(1 for line in open(file))          

print "Parameters from Simulation and PEDE:"
print "moduleN ", moduleN
print "PEDE Trials ",lineN  
print "Number of Global Par" , parN
print "Number of fixed modules", constN
print "Number of variable modules", alignN
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
		for i_par in range(0, moduleN*parN):
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

# print data 
# print trackN	

##################PLOTING##############################

plt.rcParams.update({'font.size': 14})
#Plot difference for all modules
plt.figure(1)

#Loop over expected parameters
for i_par in range(0, len(expectPars)):

	splitLabel = [int(x) for x in str(expectPars[i_par])]  # 0 = Module, 1= Parameter

	#put a line at  0 on a plot
	axes = plt.gca()
	axes.locator_params(nbins=4, axis='y')
	line = [[0,0], [trackN[lineN-1], 0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')
	#Loop over data lines 

	for i_line in range(0, lineN):
		
		dM=(data[i_par][i_line][1]-mis_C[i_par])*1e3  # mm to um rad to mrad 
		errorM=data[i_par][i_line][2]*1e3
		plt.errorbar(trackN[i_line], dM, yerr=errorM, color="red") # converting 1 cm = 10'000 um
		plt.plot(trackN[i_line], dM, marker="_", color="red")
		axes.set_xlim(-500,trackN[lineN-1]+100)
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
			plt.title('FoM M%s X'  %(int(splitLabel[0])) , fontsize=16)
			
		if(splitLabel[1]==2):
			plt.title('FoM M%s Y'  %(int(splitLabel[0])) , fontsize=16)
		
		if(splitLabel[1]==3):
			plt.title('FoM M%s $\Phi$'  %(int(splitLabel[0])) , fontsize=16)
		
		if(splitLabel[1]==4):
			plt.title('FoM M%s $\Theta$'   %(int(splitLabel[0])) , fontsize=16)

		if(splitLabel[1]==5):
			plt.title('FoM M%s $\Psi$'   %(int(splitLabel[0])) , fontsize=16)
		

	plt.savefig(str(expectPars[i_par])+".png")
	plt.clf()
		

print "Plots saved from:" , str(file) , "on", strftime("%Y-%m-%d %H:%M:%S")