#!/usr/bin/python
import itertools
import numpy as np  # smart arrays 
import subprocess
import argparse, sys

#Create combinations
combinations = np.array([''.join(x) for x in itertools.combinations('12345678',2)])

print "Total combinations:", len(combinations)


firstModuleX=[]
firstModuleY=[]
secondModuleX=[]
secondModuleY=[]
for i in range(0, len(combinations)):
	firstModuleX.append(str(combinations[i])[:1]+"1")
	firstModuleY.append(str(combinations[i])[:1]+"2")
	secondModuleX.append(str(combinations[i])[1:2]+"1")
	secondModuleY.append(str(combinations[i])[1:2]+"2")

firstModuleX=np.array(firstModuleX)
firstModuleY=np.array(firstModuleY)
secondModuleX=np.array(secondModuleX)
secondModuleY=np.array(secondModuleY)

# print firstModule
# print secondModule


#Now Loop over the combinations and produce FoM
for i in range (0, len(combinations)):
# for i in range (0, 8):

	#Write new steering file
	f = open('ParameterFile.txt', "w")
	f.write("PARAMETERS\n")
	f.write(str(firstModuleX[i]) + " 0.0 " + "-1\n")
	f.write(str(firstModuleY[i]) + " 0.0 " + "-1\n")
	f.write(str(secondModuleX[i]) + " 0.0 " + "-1\n")
	f.write(str(secondModuleY[i]) + " 0.0 " + "-1\n")
	f.write("\n")
	f.close()  


	#Run PEDE
	subprocess.call(["/Users/gleb/software/alignTrack/PEDE/pede", "SteeringFile.txt"])

	#Read new alignments
	subprocess.call(["python" , "ConcatenatePEDE.py", "-m", "w"])

	#Produce FoM
	label = "Fixed: " + str(firstModuleX[i]) + " " + str(firstModuleY[i]) + " " + str(secondModuleX[i]) + " " + str(secondModuleY[i])
	subprocess.call(["python" , "../TrackerLaunch.py", "-m", "PEDE_Mis_art.txt", "-eL", str(label)])

	#Keep a copy of the file
	subprocess.call(["cp" , "PEDE_Mis_art.txt", str(label)+"PEDE_Mis_art.txt"])
	subprocess.call(["cp" , "ParameterFile.txt", str(label)+"ParameterFile.txt"])


