####################################################################
# Appends the created list of offsets for both stations to a single
# FHICL file (for Tracking or MC) 
#
# Created: 5 November 2018 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################

import argparse, sys
import re
from operator import add
import subprocess


def getOffsets(f, name):
	copy = False #flag for processing
	offsets = [] #tmp storage buffer 
	for line in f:
	    if line.strip() == str(name)+": [": #start reading offsets after the line 
	        copy = True #found starting point 

	    elif line.strip() == "]":
	        copy = False #found ending point 

	    elif copy==True:
	        offsets.append(line.strip()) #while in between start-end add to the buffer 
	#check for par number 
	if len(offsets) != parN:
		print("Did not receive the expected number of parameters - check the inputs files!")
		sys.exit()
	else:
		print("Offsets stored")
		return offsets



parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-out", "--outFile", help='Output file')
parser.add_argument("-s12", "--station12", help='station 12')
parser.add_argument("-s18", "--station18", help='station 18')
args = parser.parse_args()
outFile = str(args.outFile)
s12File = str(args.station12) 
s18File = str(args.station18)

#Tracker constants
moduleN=8
strawN=128
parN=moduleN*strawN

# Open the output file either MC or Tracking
# and check that no previous offsets have been written,
with open(outFile) as f:
     if "services.Geometry.strawtracker.radialOffset" in f.read():
         print("The output file already contains offsets - check the correct file is passed/manually backup and delete old offsets!")
         sys.exit()
     if "services.Geometry.strawtracker.verticalOffset" in f.read():
         print("The output file already contains offsets - check the correct file is passed/manually backup and delete old offsets!")
         sys.exit()

#Store the radial and vertical offsets from S12 and S18 as lists 
s12Rad=()
s12Ver=()
s18Rad=()
s18Ver=()

f=open(s12File, "r")
s12Rad=getOffsets(f, "radialOffset")
f=open(s12File, "r")
s12Ver=getOffsets(f, "verticalOffset")
f=open(s18File, "r")
s18Rad=getOffsets(f, "radialOffset")
f=open(s18File, "r")
s18Ver=getOffsets(f, "verticalOffset")

#a+ allows to append at the end of the file
f = open(outFile, 'a+')

f.write("\n\n//Straw Offsets\n")

f.write("services.Geometry.strawtracker.radialOffset: [\n")
#for all but last entry (have commas already)
for offset in s12Rad[:-1]:
	f.write(offset+"\n")
f.write(s12Rad[-1]+",\n")   # need a ',' here to join rad and ver in one FHICL array 
for offset in s18Rad:
	f.write(offset+"\n")   #the last element doesn't need a comma 

f.write("]")

f.write("\n\n") 

f.write("services.Geometry.strawtracker.verticalOffset: [\n")
for offset in s12Ver[:-1]:
	f.write(offset+"\n")
f.write(s12Ver[-1]+",\n")     # need a ',' here to join rad and ver in one FHICL array 
for offset in s18Ver:
	f.write(offset+"\n") #the last element doesn't need a comma 
	
f.write("]")

f.close()