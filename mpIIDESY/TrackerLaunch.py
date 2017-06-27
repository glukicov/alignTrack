#!/usr/bin/python

####################################################################
# FoM Plots for comparison of actual misalignment vs PEDE results 
#
# 
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################


import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines 


# Getting constants from MC
with open("Tracker_p_constants.txt") as f:
     for line in f:  #Line is a string
        number_str = line.split()
        moduleN=int(number_str[0])

# Quickly open the PEDe file and count lines only:
lineN= sum(1 for line in open('PEDE_Mis.txt'))
            

print "Parameters from Simulation and PEDE:"
print "moduleN= ",moduleN
print "PEDE Trials= ",lineN  


mis_C = [0 for i_module in xrange(moduleN)]
# Get 1 set of misalignment from simulation 

with open("Tracker_pede_mis.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
        
        for i_module in range(0, moduleN):
        	mis_C[i_module]=abs(float(number_str[i_module]))
     

Labels = [[0 for i_module in xrange(moduleN)] for i_lines in xrange(lineN)] 
Misals = [[0 for i_module in xrange(moduleN)] for i_lines in xrange(lineN)] 
Errors = [[0 for i_module in xrange(moduleN)] for i_lines in xrange(lineN)] 
trackN = [] # track count correspond to line number 


with open("PEDE_Mis.txt") as f:
	line_i = 0
	for line in f:  #Line is a string
		number_str = line.split()
		
		for i_module in range(0, moduleN):
			Labels[line_i][i_module]=int(number_str[0+i_module*3])
			Misals[line_i][i_module]=abs(float(number_str[1+i_module*3]))
			Errors[line_i][i_module]=float(number_str[2+i_module*3])
	
		trackN.append(int(number_str[12]))
		line_i = line_i + 1
		


##################PLOTING##############################


# for i_module in range(0, moduleN):
# 	for i_lines in range(0, lineN):
# 		dM=Misals[i_lines][i_module] - mis_C[i_module]
# 		print dM

# '''

#Plot un-normalised difference for all modules
# plt.figure(1)
# for i_module in range(0, moduleN):
		
# 	if(i_module==0):
# 		plotID=221
# 	if(i_module==1):
# 		plotID=222
# 	if(i_module==2):
# 		plotID=223
# 	if(i_module==3):
# 		plotID=224
# 	if(i_module==4):
# 		plotID=231
# 	if(i_module==5):
# 		plotID=213

# 	plt.subplot(plotID)
# 	axes = plt.gca()
# 	line = [[0,0], [trackN[lineN-1], 0]]
# 	plt.plot(
# 	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
# 	    color = 'green')
# 	for i_lines in range(0, lineN):
# 		#dM=Misals[i_lines][i_module] - mis_C[i_module]
# 		dM=Misals[i_lines][i_module]-mis_C[i_module]
# 		#print i_lines, i_module, dM
# 		plt.errorbar(trackN[i_lines], dM*1e4, yerr=Errors[i_lines][i_module]*1e4, color="red") # converting 1 cm = 10'000 um
# 		plt.xlabel("Number of Tracks")
# 		plt.ylabel("$\delta$ Misalignment [um]")
# 		plt.title('FoM Module %s' %(i_module))
# 		#axes.set_ylim([beamX0-1,beamX1+1])
# 		axes.set_xlim(-1000,trackN[lineN-1])

# plt.show()


comparatorError=[]
comparatorM=[]

Sum_Chi2_1=0
Sum_Chi2_2=0

plt.figure(1)
for i_module in range(1, moduleN-1):
			
	if(i_module==1):
		plotID=211
	if(i_module==2):
		plotID=212
	if(i_module==3):
		plotID=224
	if(i_module==4):
		plotID=231

	plt.subplot(plotID)
	axes = plt.gca()
	line = [[0,0], [trackN[lineN-1], 0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')
	for i_lines in range(0, lineN):
		dM=(Misals[i_lines][i_module]-mis_C[i_module])*1e4
		errorM=Errors[i_lines][i_module]*1e4
		print "i_module=", i_module, "dM= ", dM, "error= ", errorM
		comparatorM.append(abs(dM))
		comparatorError.append(abs(errorM))
		plt.errorbar(trackN[i_lines], dM, yerr=errorM, color="red") # converting 1 cm = 10'000 um
		plt.xlabel("Number of Tracks")
		plt.ylabel("$\Delta$ Misalignment [um]")
		plt.title('FoM Module %s' %(i_module))
		#axes.set_ylim([beamX0-1,beamX1+1])
		axes.set_xlim(-1000,trackN[lineN-1])
		smalldM=min(comparatorM)
		hugedM=max(comparatorM)
		hugeError=max(comparatorError)
		smallError=min(comparatorError)
		
	
	print "Smallest Error=", smallError, "um in Module", i_module
	print "Largest Error=", hugeError, "um in Module", i_module
	print "Smallest dM=", smalldM, "um in Module", i_module
	print "Largest dM=", hugedM, "um in Module", i_module
	del comparatorM[:]
	del comparatorError[:]

plt.show()

Chi2ndf_1=[]
Chi2ndf_2=[]
Tracks1=[]
Tracks2=[]

for i_module in range(1, moduleN-1):
	for i_cut in range(0, lineN-4):
		for i_lines in range(i_cut, lineN):
			dM=(Misals[i_lines][i_module]-mis_C[i_module])*1e4
			errorM=Errors[i_lines][i_module]*1e4
			if (i_module==1):
				Sum_Chi2_1=Sum_Chi2_1+(dM*dM)/(errorM*errorM)
				#print Sum_Chi2_1
			if (i_module==2):	
				Sum_Chi2_2=Sum_Chi2_2+(dM*dM)/(errorM*errorM)
		ndf=lineN-1-i_cut
		
		#print "ndf=", ndf, "at cut ", Tracks[i_cut]
		if (i_module==1):
			Chi2ndf_1.append(Sum_Chi2_1/ndf)
			Tracks1.append(float(trackN[i_cut]))
		if (i_module==2):
			Chi2ndf_2.append(Sum_Chi2_2/ndf)
			Tracks2.append(float(trackN[i_cut]))
		
		Sum_Chi2_1=0
		Sum_Chi2_2=0
		


plt.figure(2)			
if(i_module==1):
	plotID=211
if(i_module==2):
	plotID=212
plt.subplot(plotID)
axes = plt.gca()
for i_items in range(0, len(Tracks1)-1):
	plt.plot(Tracks1[i_items], Chi2ndf_1[i_items], color="red", marker="x")
	#plt.plot(Tracks2[i_items] ,Chi2ndf_2[i_items], color="green", marker="x")

plt.xlabel("Exclusive Track-cut (>= # Tracks)")
plt.ylabel("Chi2 [um]")
plt.title("Chi^2 for M1 (red) and M2 (green)")
axes.set_xlim(-1000,trackN[lineN-1])
plt.show()





		
