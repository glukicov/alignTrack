#!/usr/bin/python

####################################################################
# Sanity plots for Tracker Alignment.  
# FoM Plots for comparison of actual misalignment vs PEDE results 
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 8 January 2018 by Gleb
#####################################################################
import argparse, sys

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--mode', help='mode')
args = parser.parse_args()

file = str(args.mode)

import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines 

# # Getting constants from MC
# with open("Tracker_p_constants.txt") as f:
#      for line in f:  #Line is a string
#         number_str = line.split()
#         moduleN=int(number_str[0])

parN=5 # GL
constN = 2 # TODO 
moduleN = 4 
# expectPars=(21, 31)
# expectPars=(21, 22, 31, 32)
#expectPars=(21, 22, 23, 31, 32, 33)
expectPars=(21, 22, 23, 24, 25, 31, 32, 33, 34, 35)

# Quickly open the PEDe file and count lines only:
lineN= sum(1 for line in open(file))
            

print "Parameters from Simulation and PEDE:"
print "moduleN= ", moduleN
print "PEDE Trials= ",lineN  
print "Number of Global Par=" , parN
print "Number of fixed modules", constN

# mis_C = (0, -0.0122, 0.00873, 0) # Theta
# mis_C = (0, 0.0087, -0.0070, 0) # PSi 
# mis_C = (0.0, 0.2, -0.1, 0.0) # Y 
# mis_C = (0.0, 0.0, 0.15, 0.2, -0.2, -0.1, 0.0, 0.0)  # XY
# mis_C = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0,  -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,)
# mis_C = (0.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.0, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)  // X
# mis_C = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0087, 0.0, 0.0, 0.0, 0.0, -0.0175, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) // Phi 
# mis_C = (0.0, 0.0, 0.0, 0.00872665, 0.00872665, 0.00872665, -0.00872665, -0.00872665, -0.00872665, 0.0, 0.0, 0.0)
# mis_C = (0.0, 0.0, 0.0, 0.0, 0.0, 0.2, -0.1, 0.0087, -0.0122, 0.0087, -0.1, 0.15,  -0.0175, 0.00873, -0.0070, 0.0, 0.0, 0.0, 0.0, 0.0)

# mis_C = ( 0.0, 0.0, 0.0, 0.0, 0.0,   0.2, -0.1, 0.0087, -0.0122, 0.0087,    -0.1, 0.15, -0.0175, 0.00873, -0.0070,   0.0, 0.0, 0.0, 0.0, 0.0)
mis_C = ( 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0)

print "Misal"

# mis_C = [0 for i_module in xrange((moduleN-constN)*parN)]
# Get 1 set of misalignment from simulation 

# with open("Tracker_pede_mis.txt") as f:
# 	for line in f:  #Line is a string
# 		number_str = line.split()
        
#         for i_module in range(0, (moduleN-constN)*parN):
#         	mis_C[i_module]=float(number_str[i_module])


Labels = [[0 for i_module in xrange((moduleN)*parN)] for i_lines in xrange(lineN)] 
Misals = [[0 for i_module in xrange((moduleN)*parN)] for i_lines in xrange(lineN)] 
Errors = [[0 for i_module in xrange((moduleN)*parN)] for i_lines in xrange(lineN)] 
trackN = [] # track count correspond to line number 


with open(file) as f:
	line_i = 0
	for line in f:  #Line is a string
		number_str = line.split()

		for i_par in range(0, moduleN*parN):
			label=int(number_str[0+i_par*3])
			#if (label==21 or label==22 or label==31 or label==32):
			# if (label==21 or label==22 or  label==23 or label==31 or label==32 or label==33):
			if (label==21 or label==22 or  label==23 or  label==24 or label==25 or label==31 or label==32 or label==33 or label==34 or label==35):
			# if (label==21 or label==31):
				misal=float(number_str[1+i_par*3])
				error=float(number_str[2+i_par*3])
				# print misal
				# print error
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
		

# with open("Tracker_metric.txt") as f:
# 	line_i = 0
# 	for line in f:  #Line is a string
# 		metric = line

print "PEDE Mis: ", Misals
print "PEDE Errors:" , Errors
print "Truth Misalignment: " ,mis_C

##################PLOTING##############################


# for i_module in range(0, moduleN):
# 	for i_lines in range(0, lineN):
# 		dM=Misals[i_lines][i_module] - mis_C[i_module]
# 		print dM

plt.rcParams.update({'font.size': 8})

moduleN=moduleN-constN  # Only look at misaligned modules 

#Plot difference for all modules
plt.figure(1)
plots = []
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
			plotID=223
		if(i_par==31):
			plotID=222
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

	if (parN==5):
		if(i_par==21):
			plotID=221
		if(i_par==22):
			plotID=223
		# if(i_par==23):
		# 	plotID=525
		# if(i_par==24):
		# 	plotID=527
		# if(i_par==25):
		# 	plotID=527
		# if(i_par==31):
		# 	plotID=522
		# if(i_par==32):
		# 	plotID=524
		# if(i_par==33):
		# 	plotID=526
		# if(i_par==34):
		# 	plotID=528
		# if(i_par==35):
		# 	plotID=528

	plt.subplot(plotID)
	plt.tight_layout()
	axes = plt.gca()
	axes.locator_params(nbins=4, axis='y')
	line = [[0,0], [trackN[lineN-1], 0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')
	for i_lines in range(0, lineN):
		if(i_par==21 or i_par==31):
		# if(i_par==21 or i_par==31 or i_par==22 or i_par==32 or i_par==23 or i_par==33 or i_par==24 or i_par==34):
		#if(i_par==21 or i_par==31 or i_par==22 or i_par==32 or i_par==23 or i_par==33 or i_par==25 or i_par==35):
			dM=(Misals[i_lines][i_counter+parN]-mis_C[i_counter+parN])*1e3  # mm to um rad to mrad 
			#dM=(Misals[i_lines][i_counter+2]-mis_C[i_counter])*1e4
			#print 'Misals[i_lines][i_module]=', Misals[i_lines][i_module], 'mis_C[i_module]=', mis_C[i_module], 'dM=', dM
			errorM=Errors[i_lines][i_counter+parN]*1e3
			#errorM=Errors[i_lines][i_counter+2]*1e4
			print "dM=", dM, "errorM=", errorM
			plt.errorbar(trackN[i_lines], dM, yerr=errorM, color="red") # converting 1 cm = 10'000 um
			plt.plot(trackN[i_lines], dM, marker="_", color="red")
			
			#axes.set_ylim([beamX0-1,beamX1+1])
			axes.set_xlim(-500,trackN[lineN-1]+100)
		
		if(i_par==21 or i_par==31):
			axes.set_ylim(-150, 150)
			plt.title('FoM M%s X'  %(int(i_par/10)) , fontsize=6)
			plt.ylabel("$\Delta$ Misalignment [um]", fontsize=6)

		# if(i_par==22 or i_par==32):
		# 	axes.set_ylim(-150, 150)
		# 	plt.title('FoM M%s Y'  %(int(i_par/10)) , fontsize=6)
		# 	plt.ylabel("$\Delta$ Misalignment [um]", fontsize=6)

		# if(i_par==23 or i_par==33):
		# 	axes.set_ylim(-10, 10)
		# 	plt.title('FoM M%s $\Phi$'  %(int(i_par/10)) , fontsize=6)
		# 	plt.ylabel("$\Delta$ Misalignment [mrad]", fontsize=6)

		# if(i_par==24 or i_par==34):
		# 	axes.set_ylim(-10, 10)
		# 	plt.title('FoM M%s $\Theta$'  %(int(i_par/10)) , fontsize=6)
		# 	plt.ylabel("$\Delta$ Misalignment [mrad]", fontsize=6)

		# if(i_par==25 or i_par==35):
		# 	axes.set_ylim(-2, 2)
		# 	plt.title('FoM M%s $\Psi$'  %(int(i_par/10)) , fontsize=6)
		# 	plt.ylabel("$\Delta$ Misalignment [mrad]", fontsize=6)

		plt.xlabel("Number of Tracks", fontsize=5)
plt.subplots_adjust(hspace=.7)
plt.gcf().subplots_adjust(top=0.90)
#plt.suptitle(str(metric), fontsize=6, style='oblique',  color="green")
plt.savefig(str(file)+".png")
print "File produced:" , str(file), ".png" 