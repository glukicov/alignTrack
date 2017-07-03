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
	
		if (moduleN==4):
			trackN.append(int(number_str[12]))
		if (moduleN==6):
		    trackN.append(int(number_str[18]))
				#print number_str[18]
		line_i = line_i + 1
		


##################PLOTING##############################


# for i_module in range(0, moduleN):
# 	for i_lines in range(0, lineN):
# 		dM=Misals[i_lines][i_module] - mis_C[i_module]
# 		print dM

if (moduleN==6):
	plt.rcParams.update({'font.size': 8})

#Plot difference for all modules
plt.figure(1)
for i_module in range(0, moduleN):
		
	if (moduleN==4):
		if(i_module==0):
			plotID=221
		if(i_module==1):
			plotID=222
		if(i_module==2):
			plotID=223
		if(i_module==3):
			plotID=224
	if (moduleN==6):
		if(i_module==0):
			plotID=421
		if(i_module==1):
			plotID=422
		if(i_module==2):
			plotID=423
		if(i_module==3):
			plotID=424
		if(i_module==4):
			plotID=425
		if(i_module==5):
			plotID=426

	plt.subplot(plotID)
	plt.tight_layout()
	axes = plt.gca()
	axes.locator_params(nbins=4, axis='y')
	line = [[0,0], [trackN[lineN-1], 0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')
	for i_lines in range(0, lineN):
		dM=(Misals[i_lines][i_module]-mis_C[i_module])*1e4
		errorM=Errors[i_lines][i_module]*1e4
		plt.errorbar(trackN[i_lines], dM, yerr=errorM, color="red") # converting 1 cm = 10'000 um
		plt.plot(trackN[i_lines], dM, marker="_", color="red")
		
		
		plt.title('FoM Module %s' %(i_module))
		#axes.set_ylim([beamX0-1,beamX1+1])
		axes.set_xlim(-500,trackN[lineN-1]+100)
		axes.set_ylim(-20, 40)


		plt.xlabel("Number of Tracks")
		if (i_module!=2 or i_module!=3):
			plt.ylabel("$\Delta$ Misalignment [um]")
plt.subplots_adjust(hspace=.9)
plt.gcf().subplots_adjust(bottom=0.15)

plt.savefig("FoM_All.png")
#plt.show()

comparatorError=[]
comparatorM=[]

Sum_Chi2_1=0
Sum_Chi2_2=0

plt.figure(1)
for i_module in range(1, moduleN-1):
			
	
	if (moduleN==4):
		if(i_module==1):
			plotID=211
		if(i_module==2):
	 		plotID=212
	if (moduleN==6):	
		if(i_module==1):
			plotID=221
		if(i_module==2):
			plotID=222
		if(i_module==3):
			plotID=223
		if(i_module==4):
			plotID=224
		

	plt.subplot(plotID)
	plt.tight_layout()
	axes = plt.gca()
	line = [[0,0], [trackN[lineN-1], 0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')
	for i_lines in range(0, lineN):
		dM=(Misals[i_lines][i_module]-mis_C[i_module])*1e4
		errorM=Errors[i_lines][i_module]*1e4
		#print "i_module=", i_module, "dM= ", dM, "error= ", errorM
		comparatorM.append(abs(dM))
		comparatorError.append(abs(errorM))
		plt.errorbar(trackN[i_lines], dM, yerr=errorM, color="red") # converting 1 cm = 10'000 um
		plt.plot(trackN[i_lines], dM,  marker="_", color="red")
		
		plt.title('FoM Module %s' %(i_module))
		#axes.set_ylim([beamX0-1,beamX1+1])
		axes.set_xlim(-500,trackN[lineN-1]+100)
		axes.set_ylim(-20, 40)
		smalldM=min(comparatorM)
		hugedM=max(comparatorM)
		hugeError=max(comparatorError)
		smallError=min(comparatorError)
		
	
	print "In Module", i_module,": Smallest Error=", smallError, "um, Largest Error=", hugeError, "um, Smallest dM=", smalldM, "um Largest dM=", hugedM, "um"
	del comparatorM[:]
	del comparatorError[:]
	
	if (i_module!=2 or i_module!=3):
		plt.ylabel("$\Delta$ Misalignment [um]")
	
	plt.xlabel("Number of Tracks")
plt.subplots_adjust(hspace=.5)
plt.gcf().subplots_adjust(bottom=0.15)
#plt.show()
plt.savefig("FoM_Mis.png")

print("Files produced: FoM_Mis.png and FoM_All.png")




		
