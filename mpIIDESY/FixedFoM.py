#!/usr/bin/python
import itertools
import numpy as np  # smart arrays 
import subprocess
import argparse, sys

label = "Mean effect of Fixed M2: 17; 18; 27; 28; M3:135, 137, 156, 158, 167, 168, 257, 257, 268"
subprocess.call(["mv" , "PEDE_Mis_art.txt", "BK_PEDE_Mis_art.txt"])

#####2 Modules######
firstModuleX=[11, 11, 21, 21]
firstModuleY=[12, 12, 22, 22]
secondModuleX=[71, 81, 71, 81]
secondModuleY=[72, 82, 72, 82]
fileLabel="M2"
################

firstModuleX=np.array(firstModuleX)
firstModuleY=np.array(firstModuleY)

#Now Loop over the combinations and produce FoM
for i in range (0, len(firstModuleX)):

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

	#Read new alignments and append
	subprocess.call(["python" , "../ConcatenatePEDE.py", "-m", "a"])

	#Keep a copy of the files
	subprocess.call(["cp" , "ParameterFile.txt", str(fileLabel)+str(i)+"ParameterFile.txt"])
	subprocess.call(["cp" , "millepede.res", str(fileLabel)+str(i)+"millepede.res"])

# #####3 Modules######
# firstModuleX=  [11, 11, 11, 11, 11, 11, 21, 21, 21]
# firstModuleY=  [12, 12, 12, 12, 12, 12, 22, 22, 22]
# secondModuleX= [31, 31, 51, 51, 61, 61, 51, 61, 61]
# secondModuleY= [32, 32, 52, 52, 62, 62, 52, 62, 62]
# thirdModuleX=  [51, 71, 61, 81, 71, 81, 71, 71, 81]
# thirdModuleY=  [52, 72, 62, 82, 72, 82, 72, 72, 82]
# fileLabel="M3"
# ################

# firstModuleX=np.array(firstModuleX)
# firstModuleY=np.array(firstModuleY)
# secondModuleX=np.array(secondModuleX)
# secondModuleY=np.array(secondModuleY)
# thirdModuleX=np.array(thirdModuleX)
# thirdModuleY=np.array(thirdModuleY)

# #Now Loop over the combinations and produce FoM
# for i in range (0, len(firstModuleX)):

# 	#Write new steering file
# 	f = open('ParameterFile.txt', "w")
# 	f.write("PARAMETERS\n")
# 	f.write(str(firstModuleX[i]) + " 0.0 " + "-1\n")
# 	f.write(str(firstModuleY[i]) + " 0.0 " + "-1\n")
# 	f.write(str(secondModuleX[i]) + " 0.0 " + "-1\n")
# 	f.write(str(secondModuleY[i]) + " 0.0 " + "-1\n")
# 	f.write(str(thirdModuleX[i]) + " 0.0 " + "-1\n")
# 	f.write(str(thirdModuleY[i]) + " 0.0 " + "-1\n")
# 	f.write("\n")
# 	f.close()  

# 	#Run PEDE
# 	subprocess.call(["/Users/gleb/software/alignTrack/PEDE/pede", "SteeringFile.txt"])

# 	#Read new alignments and append
# 	subprocess.call(["python" , "ConcatenatePEDE.py", "-m", "a"])

# 	#Keep a copy of the files
# 	subprocess.call(["cp" , "ParameterFile.txt", str(fileLabel)+str(i)+"ParameterFile.txt"])
# 	subprocess.call(["cp" , "millepede.res", str(fileLabel)+str(i)+"millepede.res"])

# #####4 Modules######
# firstModuleX=[11, 11, 31, 21]
# firstModuleY=[12, 12, 32, 22]
# secondModuleX=[71, 61, 61, 61]
# secondModuleY=[72, 62, 62, 62]
# thirdModuleX=[81, 71, 81, 71]
# thirdModuleY=[82, 72, 82, 72]
# label = "Mean effect of Fixed M: 1-2-7-8, 1-6-7, 3-6-8, 2-6-7"
# ################




#Produce FoM
subprocess.call(["python" , "../FixedTrackerLaunch.py", "-m", "PEDE_Mis_art.txt", "-eL", str(label)])