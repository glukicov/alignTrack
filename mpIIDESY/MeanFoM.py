#!/usr/bin/python
import itertools
import numpy as np  # smart arrays 
import subprocess
import argparse, sys

#Create combinations
# combinations = np.array([''.join(x) for x in itertools.combinations('12345678',2)])

# print "Total combinations:", len(combinations)


firstModuleX=[11, 11, 21, 21]
firstModuleY=[12, 12, 22, 22]
secondModuleX=[71, 81, 71, 81]
secondModuleY=[72, 82, 72, 82]

# firstModuleXOffsets=[-0.2, -0.2, 0.08, 0.08]
# firstModuleYOffsets=[0.1, 0.1, 0.15, 0.15]
# secondModuleXOffsets=[0.2, -0.06, 0.2, -0.06]
# secondModuleYOffsets=[0.07, 0.06, 0.07, 0.06]

firstModuleXOffsets=[-0.2, -0.2, 0.08, 0.08]
firstModuleYOffsets=[0.1, 0.1, 0.15, 0.15]
secondModuleXOffsets=[0.2, -0.06, 0.2, -0.06]
secondModuleYOffsets=[0.07, 0.06, 0.07, 0.06]


label = "Mean effect of Fixed M: 1-7, 1-8, 2-7, 2-8"

# for i in range(0, len(combinations)):
# 	firstModuleX.append(str(combinations[i])[:1]+"1")
# 	firstModuleY.append(str(combinations[i])[:1]+"2")
# 	secondModuleX.append(str(combinations[i])[1:2]+"1")
# 	secondModuleY.append(str(combinations[i])[1:2]+"2")

firstModuleX=np.array(firstModuleX)
firstModuleY=np.array(firstModuleY)
secondModuleX=np.array(secondModuleX)
secondModuleY=np.array(secondModuleY)

# print firstModule
# print secondModule

subprocess.call(["mv" , "PEDE_Mis_art.txt", "BK_PEDE_Mis_art.txt"])

#Now Loop over the combinations and produce FoM
for i in range (0, len(firstModuleX)):
# for i in range (0, 8):

	#Write new steering file
	f = open('ParameterFile.txt', "w")
	f.write("PARAMETERS\n")
	f.write(str(firstModuleX[i]) + " "  +str(firstModuleXOffsets[i]) +  " -1\n")
	f.write(str(firstModuleY[i]) + " "  +str(firstModuleYOffsets[i]) + " -1\n")
	f.write(str(secondModuleX[i]) + " " +str(secondModuleXOffsets[i]) +" -1\n")
	f.write(str(secondModuleY[i]) + " " +str(secondModuleYOffsets[i]) + " -1\n")
	f.write("\n")
	f.close()  

	#Run PEDE
	subprocess.call(["/Users/gleb/software/alignTrack/PEDE/pede", "SteeringFile.txt"])

	#Read new alignments
	subprocess.call(["python" , "ConcatenatePEDE.py", "-m", "a"])

	#Keep a copy of the files
	subprocess.call(["cp" , "PEDE_Mis_art.txt", str(label)+"PEDE_Mis_art.txt"])
	subprocess.call(["cp" , "ParameterFile.txt", str(label)+"ParameterFile.txt"])
	subprocess.call(["cp" , "millepede.res", str(label)+"millepede.res"])


#Produce FoM
subprocess.call(["python" , "../MeanTrackerLaunch.py", "-m", "PEDE_Mis_art.txt", "-eL", str(label)])

