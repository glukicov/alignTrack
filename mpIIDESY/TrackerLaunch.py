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
        	mis_C[i_module]=float(number_str[i_module])
     

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
			Misals[line_i][i_module]=float(number_str[1+i_module*3])
			Errors[line_i][i_module]=float(number_str[2+i_module*3])
	
		trackN.append(int(number_str[12]))
		line_i = line_i + 1
		


##################PLOTING##############################


# for i_module in range(0, moduleN):
# 	for i_lines in range(0, lineN):
# 		dM=Misals[i_lines][i_module] - mis_C[i_module]
# 		print dM

# '''


plt.figure(1)
for i_module in range(0, moduleN):
		
	if(i_module==0):
		plotID=221
	if(i_module==1):
		plotID=222
	if(i_module==2):
		plotID=223
	if(i_module==3):
		plotID=224
	if(i_module==4):
		plotID=431
	if(i_module==5):
		plotID=413
		
	plt.subplot(plotID)
	axes = plt.gca()
	for i_lines in range(0, lineN):
		#dM=Misals[i_lines][i_module] - mis_C[i_module]
		dM=abs(Misals[i_lines][i_module]-mis_C[i_module])
		plt.errorbar(trackN[i_lines], dM, yerr=Errors[i_lines][i_module], color="red")
		plt.xlabel("Number of Tracks")
		plt.ylabel(r"\delta\ Misalignment")
		plt.title('FoM Module %s' %(i_module))

#axes.set_ylim([beamX0-1,beamX1+1])
#axes.set_xlim([beamZ0-1,beamZ1+1])




plt.show()
















		
