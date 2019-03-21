####################################################################
# Run 4 runs of PEDE with M1/M2 - M7/M8 fixed at (0.0, 0.0) 
# Need to be run from a dir with Data.bin and SteeringFile.txt 
# from alignment 
#
# python3 FixedFoM.py -s 12 
#
# or -s 18 for Station 18. FixedTrackerLaunch.py will then
# produce the plots. 
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 26 November 2018 by Gleb
#####################################################################

import itertools
import numpy as np  # smart arrays 
import subprocess
import argparse, sys

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-s', '--stationN', help='station number')
parser.add_argument('-moduleN', '--moduleN', default=8)
args = parser.parse_args()

stationN = str(args.stationN)
moduleN = int(args.moduleN) 

subprocess.call(["mv" , "PEDE_Mis_art.txt", "BK_PEDE_Mis_art.txt"])

#####2 Modules fixing######

if (stationN == "12"):
	firstModuleX= [1211, 1211, 1211, 1211, 1211, 1211, 1211]
	firstModuleY= [1212, 1212, 1212, 1212, 1212, 1212, 1212]
	secondModuleX=[1221, 1231, 1241, 1251, 1261, 1271, 1281]
	secondModuleY=[1222, 1232, 1242, 1252, 1262, 1272, 1282]

if (stationN == "18"):
	firstModuleX= [1811, 1811, 1811, 1811, 1811, 1811, 1811]
	firstModuleY= [1812, 1812, 1812, 1812, 1812, 1812, 1812]
	secondModuleX=[1821, 1831, 1841, 1851, 1861, 1871, 1881]
	secondModuleY=[1822, 1832, 1842, 1852, 1862, 1872, 1882]


fileLabel="2MFixed"
################

#Fix these 2M at 0.0 
FixedPositionX1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
FixedPositionY1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
FixedPositionX2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
FixedPositionY2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# FixedPositionX1 = [-0.2, -0.2, 0.08, 0.08]
# FixedPositionY1 = [0.1, 0.1, 0.15, 0.15]
# FixedPositionX2 = [0.2, 0.2, -0.06, -0.06]
# FixedPositionY2 = [0.07, 0.07, 0.06, 0.06]


firstModuleX=np.array(firstModuleX)
firstModuleY=np.array(firstModuleY)

#Now Loop over the combinations and produce FoM
for i in range (0, len(firstModuleX)):

	#set new  label
	label = "S"+str(stationN)+"_"+str(i+2)

	#Write new steering file
	f = open('ParameterFile.txt', "w")
	f.write("PARAMETERS\n")
	f.write(str(firstModuleX[i]) + " " + str(FixedPositionX1[i]) + " -1\n")
	f.write(str(firstModuleY[i]) + " " + str(FixedPositionY1[i]) + " -1\n")
	f.write(str(secondModuleX[i]) +  " " + str(FixedPositionX2[i]) + " -1\n")
	f.write(str(secondModuleY[i]) +  " " + str(FixedPositionY2[i]) + " -1\n")
	f.write("\n")
	f.close()  

	#Run PEDE
	subprocess.call(["/Users/gleb/software/alignTrack/PEDE/pede", "SteeringFile.txt"])

	#Read new alignments and over
	subprocess.call(["python" , "../../ConcatenatePEDE.py", "-m", "w", "-moduleN", str(moduleN)])

	#Keep a copy of the files
	subprocess.call(["cp" , "ParameterFile.txt", str(fileLabel)+str(i+1)+"ParameterFile.txt"])
	subprocess.call(["cp" , "millepede.res", str(fileLabel)+str(i+1)+"millepede.res"])
	
	#Produce FoM
	subprocess.call(["python3", "../../FixedTrackerLaunch.py", "-m", "PEDE_Mis_art.txt", "-eL", str(label), "-s", str(stationN)])

   # subprocess.call(["python3", "../../RobustTrackerLaunch.py", "-m", "PEDE_Mis_art.txt", "-eL", str(label), "-s", str(stationN), "-moduleN", str(moduleN)])