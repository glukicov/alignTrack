#!/usr/bin/python
import itertools
import numpy as np  # smart arrays 
import subprocess
import argparse, sys

#Create combinations
# combinations = np.array([''.join(x) for x in itertools.combinations('12345678',2)])

# print "Total combinations:", len(combinations)

# #####1 Module######
# firstModuleX=[11, 21, 31, 41, 51, 61, 71, 81]
# firstModuleY=[12, 22, 32, 42, 52, 62, 72, 82]
# label = "Mean effect of Fixed M: 1-8"
# ################

# #####2 Modules######
# firstModuleX=[11, 11, 21, 21]
# firstModuleY=[12, 12, 22, 22]
# secondModuleX=[71, 81, 71, 81]
# secondModuleY=[72, 82, 72, 82]
# label = "Mean effect of Fixed M: 1-7, 1-8, 2-7, 2-8"
# ################


#####3 Modules######
firstModuleX=  [11, 11, 11, 11, 11, 11, 21, 21, 21, 21]
firstModuleY=  [12, 12, 12, 12, 12, 12, 22, 22, 22, 22]
secondModuleX= [31, 31, 51, 51, 61, 61, 51, 61, 61, 71]
secondModuleY= [32, 32, 51, 52, 62, 62, 52, 62, 62, 71]
thirdModuleX=  [51, 71, 61, 81, 71, 81, 71, 71, 81, 81]
thirdModuleY=  [52, 72, 62, 82, 72, 82, 72, 72, 82, 82]
label = "Mean effect of Fixed M: 135, 137, 156, 158, 167, 168, 257, 267, 268, 278"
################

# #####4 Modules######
# firstModuleX=[11, 11, 31, 21]
# firstModuleY=[12, 12, 32, 22]
# secondModuleX=[71, 61, 61, 61]
# secondModuleY=[72, 62, 62, 62]
# thirdModuleX=[81, 71, 81, 71]
# thirdModuleY=[82, 72, 82, 72]
# label = "Mean effect of Fixed M: 1-2-7-8, 1-6-7, 3-6-8, 2-6-7"
# ################


# for i in range(0, len(combinations)):
# 	firstModuleX.append(str(combinations[i])[:1]+"1")
# 	firstModuleY.append(str(combinations[i])[:1]+"2")
# 	secondModuleX.append(str(combinations[i])[1:2]+"1")
# 	secondModuleY.append(str(combinations[i])[1:2]+"2")

firstModuleX=np.array(firstModuleX)
firstModuleY=np.array(firstModuleY)
secondModuleX=np.array(secondModuleX)
secondModuleY=np.array(secondModuleY)
thirdModuleX=np.array(thirdModuleX)
thirdModuleY=np.array(thirdModuleY)

# print firstModule
# print secondModule

subprocess.call(["mv" , "PEDE_Mis_art.txt", "BK_PEDE_Mis_art.txt"])

#Now Loop over the combinations and produce FoM
for i in range (0, len(firstModuleX)):
# for i in range (0, 8):

	#Write new steering file
	f = open('ParameterFile.txt', "w")
	f.write("PARAMETERS\n")
	f.write(str(firstModuleX[i]) + " 0.0 " + "-1\n")
	f.write(str(firstModuleY[i]) + " 0.0 " + "-1\n")
	f.write(str(secondModuleX[i]) + " 0.0 " + "-1\n")
	f.write(str(secondModuleY[i]) + " 0.0 " + "-1\n")
	f.write(str(thirdModuleX[i]) + " 0.0 " + "-1\n")
	f.write(str(thirdModuleY[i]) + " 0.0 " + "-1\n")
	f.write("\n")
	f.close()  

	#Run PEDE
	subprocess.call(["/Users/gleb/software/alignTrack/PEDE/pede", "SteeringFile.txt"])

	#Read new alignments
	subprocess.call(["python" , "ConcatenatePEDE.py", "-m", "a"])

	#Keep a copy of the files
	#subprocess.call(["cp" , "PEDE_Mis_art.txt", str(i)+"PEDE_Mis_art.txt"])
	subprocess.call(["cp" , "ParameterFile.txt", str(i)+"ParameterFile.txt"])
	subprocess.call(["cp" , "millepede.res", str(i)+"millepede.res"])

#Produce FoM
subprocess.call(["python" , "../FixedTrackerLaunch.py", "-m", "PEDE_Mis_art.txt", "-eL", str(label)])

