#!/usr/bin/python
import itertools
import numpy as np  # smart arrays 
import subprocess
import argparse, sys

#Create combinations
combinations = np.array([''.join(x) for x in itertools.combinations('12345678',2)])

print "Total combinations:", len(combinations)


firstModule=[]
secondModule=[]
for i in range(0, len(combinations)):
	firstModule.append(str(combinations[i])[:1]+"1")
	secondModule.append(str(combinations[i])[1:2]+"1")

firstModule=np.array(firstModule)
secondModule=np.array(secondModule)

# print firstModule
# print secondModule


#Now Loop over the combinations and produce FoM
for i in range (0, len(combinations)):
# for i in range (0, 8):

	#Write new steering file
	f = open('ParameterFile.txt', "w")
	f.write("PARAMETERS\n")
	f.write(str(firstModule[i]) + " 0.0 " + "-1\n")
	f.write(str(secondModule[i]) + " 0.0 " + "-1\n")
	f.write("\n")
	f.close()  


	#Run PEDE
	subprocess.call(["/Users/gleb/software/alignTrack/PEDE/pede", "SteeringFile.txt"])

	#Read new alignments
	subprocess.call(["python" , "ConcatenatePEDE.py", "-m", "w"])

	#Produce FoM
	label = "Fixed: " + str(firstModule[i]) + " " + str(secondModule[i])
	subprocess.call(["python" , "../TrackerLaunch.py", "-m", "PEDE_Mis_art.txt", "-eL", str(label)])

	#Keep a copy of the file
	subprocess.call(["cp" , "PEDE_Mis_art.txt", str(label)+"PEDE_Mis_art.txt"])
	subprocess.call(["cp" , "ParameterFile.txt", str(label)+"ParameterFile.txt"])


