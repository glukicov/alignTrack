#!/usr/bin/python

####################################################################
# FoM Plots for comparison of actual misalignment vs PEDE results 
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

parN=3
constN =2 # TODO 
expectPars=(21, 22, 23, 31, 32, 33)
# Quickly open the PEDe file and count lines only:
lineN= sum(1 for line in open('PEDE_Mis.txt'))
            

print "Parameters from Simulation and PEDE:"
print "moduleN= ", moduleN
print "PEDE Trials= ",lineN  


mis_C = [0 for i_module in xrange((moduleN-constN)*parN)]
# Get 1 set of misalignment from simulation 

with open("Tracker_pede_mis.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
        
        for i_module in range(0, (moduleN-constN)*parN):
        	mis_C[i_module]=float(number_str[i_module])



Labels = [[0 for i_module in xrange((moduleN)*parN)] for i_lines in xrange(lineN)] 
Misals = [[0 for i_module in xrange((moduleN)*parN)] for i_lines in xrange(lineN)] 
Errors = [[0 for i_module in xrange((moduleN)*parN)] for i_lines in xrange(lineN)] 
trackN = [] # track count correspond to line number 


with open("PEDE_Mis.txt") as f:
	line_i = 0
	for line in f:  #Line is a string
		number_str = line.split()

		for i_par in range(0, moduleN*parN):
			label=int(number_str[0+i_par*3])
			# if (label==21 or label==22 or label==31 or label==32):
			if (label==21 or label==22 or  label==23 or label==31 or label==32 or label==33):
				misal=float(number_str[1+i_par*3])
				error=float(number_str[2+i_par*3])
				Labels[line_i][i_par]=label
				Misals[line_i][i_par]=misal
				Errors[line_i][i_par]=error
				#print "label=", label, "misal=", misal, "error=", error
			else:
				i_par+1
	
		if (moduleN==4):
			trackN.append(int(number_str[-1]))
		if (moduleN==6):
		    trackN.append(int(number_str[-1]))
				#print number_str[18]
		line_i = line_i + 1
		

with open("Tracker_metric.txt") as f:
	line_i = 0
	for line in f:  #Line is a string
		metric = line

print Misals
print Errors
print mis_C

##################PLOTING##############################


# for i_module in range(0, moduleN):
# 	for i_lines in range(0, lineN):
# 		dM=Misals[i_lines][i_module] - mis_C[i_module]
# 		print dM

plt.rcParams.update({'font.size': 8})

moduleN=moduleN-constN  

#Plot difference for all modules
plt.figure(1)
for i_counter in range(0, parN*moduleN):
	i_par=expectPars[i_counter]
	if (parN==1):
		if(i_par==21):
			plotID=211
		if(i_par==31):
			plotID=212
		
	if (parN==2):
		if(i_par==21):
			plotID=221
		if(i_par==22):
			plotID=222
		if(i_par==31):
			plotID=223
		if(i_par==32):
			plotID=224

	if (parN==3):
		if(i_par==21):
			plotID=321
		if(i_par==22):
			plotID=323
		if(i_par==23):
			plotID=325
		if(i_par==31):
			plotID=322
		if(i_par==32):
			plotID=324
		if(i_par==33):
			plotID=326

	plt.subplot(plotID)
	plt.tight_layout()
	axes = plt.gca()
	axes.locator_params(nbins=4, axis='y')
	line = [[0,0], [trackN[lineN-1], 0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')
	for i_lines in range(0, lineN):
		dM=(Misals[i_lines][i_counter+3]-mis_C[i_counter])*1e4 # cm -> um 
		if(i_par==23 or i_par==33):
			dM=(Misals[i_lines][i_counter+3]-mis_C[i_counter])*1e3  #rad -> urad 
		#dM=(Misals[i_lines][i_counter+2]-mis_C[i_counter])*1e4
		#print 'Misals[i_lines][i_module]=', Misals[i_lines][i_module], 'mis_C[i_module]=', mis_C[i_module], 'dM=', dM
		errorM=Errors[i_lines][i_counter+3]*1e4
		#errorM=Errors[i_lines][i_counter+2]*1e4
		print "dM=", dM, "errorM=", errorM
		plt.errorbar(trackN[i_lines], dM, yerr=errorM, color="red") # converting 1 cm = 10'000 um
		plt.plot(trackN[i_lines], dM, marker="_", color="red")
		
		#axes.set_ylim([beamX0-1,beamX1+1])
		axes.set_xlim(-500,trackN[lineN-1]+100)
		if(i_par==21 or i_par==31):
			axes.set_ylim(-10, 10)
			plt.title('FoM M%s x'  %(int(i_par/10)) , fontsize=10)
		if(i_par==22 or i_par==32):
			axes.set_ylim(-2000, 2000)
			plt.title('FoM M%s z'  %(int(i_par/10)) , fontsize=10)
		if(i_par==23 or i_par==33):
			axes.set_ylim(-10, 10)
			plt.title('FoM M%s $\Theta$'  %(int(i_par/10)) , fontsize=10)

		plt.xlabel("Number of Tracks", fontsize=10)
		if (i_par==21 or i_par==31 or i_par==22 or i_par==32):
			plt.ylabel("$\Delta$ Misalignment [um]")
		if (i_par==23 or i_par==33):
			plt.ylabel("$\Delta$ Misalignment [urad]")
plt.subplots_adjust(hspace=.4)
plt.gcf().subplots_adjust(top=0.90)
plt.suptitle(str(metric), fontsize=6, style='oblique',  color="green")
plt.savefig("FoM_All.png")
print("File produced: FoM_All.png")