#!/usr/bin/python
import itertools
import numpy as np  # smart arrays 
import subprocess
import argparse, sys

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-s', '--stationN', help='station number')
args = parser.parse_args()

stationN = str(args.stationN)

label = "S"+str(stationN)
subprocess.call(["mv" , "PEDE_Mis_art.txt", "BK_PEDE_Mis_art.txt"])

#####2 Modules######
if (stationN == "10"):
	firstModuleX=[1011,  1011, 1021, 1021]
	firstModuleY=[1012,  1012, 1022, 1022]
	secondModuleX=[1071, 1081, 1071, 1081]
	secondModuleY=[1072, 1082, 1072, 1082]

if (stationN == "18"):
	firstModuleX= [1811,  1811, 1821, 1821]
	firstModuleY= [1812,  1812, 1822, 1822]
	secondModuleX=[1871,  1881, 1871, 1881]
	secondModuleY=[1872,  1882, 1872, 1882]

if (stationN == "12"):
	firstModuleX= [1211,  1211, 1221, 1221]
	firstModuleY= [1212,  1212, 1222, 1222]
	secondModuleX=[1271,  1281, 1271, 1281]
	secondModuleY=[1272,  1282, 1272, 1282]

fileLabel="2MFixed"
################

FixedPositionX1 = [0.0, 0.0, 0.0, 0.0]
FixedPositionY1 = [0.0, 0.0, 0.0, 0.0]
FixedPositionX2 = [0.0, 0.0, 0.0, 0.0]
FixedPositionY2 = [0.0, 0.0, 0.0, 0.0]

# FixedPositionX1 = [-0.2, -0.2, 0.08, 0.08]
# FixedPositionY1 = [0.1, 0.1, 0.15, 0.15]
# FixedPositionX2 = [0.2, 0.2, -0.06, -0.06]
# FixedPositionY2 = [0.07, 0.07, 0.06, 0.06]


firstModuleX=np.array(firstModuleX)
firstModuleY=np.array(firstModuleY)

#Now Loop over the combinations and produce FoM
for i in range (0, len(firstModuleX)):

	#set new 

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

	#Read new alignments and append
	subprocess.call(["python" , "../../ConcatenatePEDE.py", "-m", "a"])

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
subprocess.call(["../../FixedTrackerLaunch.py", "-m", "PEDE_Mis_art.txt", "-eL", str(label), "-s", str(stationN)])