#!/usr/bin/env python

import os
import string
import time
import decimal
from ROOT import *

Good_C_hits=[]
Bad_C_hits=[]


#Read file and store in lists 
with open("../C_mpIIDESY/Good_C_Mp2debug_mp2.txt") as f:
    #next(f)
    for lineC in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number = lineC.split()
        Good_C_hits.append(number[0])
           
        
   
#Read file and store in lists 
with open("../C_mpIIDESY/Bad_C_Mp2debug_mp2.txt") as f:
    #next(f)
    for lineC in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number = lineC.split()
        Bad_C_hits.append(number[0])
'''

with open("../C_mpIIDESY/Good_C_Mp2debug_off.txt") as f:
    #next(f)
    for lineC in f:  
        number = lineC.split()
        

        
with open("../C_mpIIDESY/Bad_C_Mp2debug_off.txt") as f:
    #next(f)
    for lineC in f:  
        number = lineC.split()
        


with open("../C_mpIIDESY/Good_C_Mp2debug_calc.txt") as f:
    #next(f)
    for lineC in f:
        number = lineC.split()
        
        

with open("../C_mpIIDESY/Bad_C_Mp2debug_calc.txt") as f:
    #next(f)
    for lineC in f:
        number = lineC.split()

with open("../C_mpIIDESY/Good_C_Mp2debug_MC.txt") as f:
    #next(f)
    for lineC in f:
        number = lineC.split()
    

with open("../C_mpIIDESY/Bad_C_Mp2debug_MC.txt") as f:
    #next(f)
    for lineC in f:
        number = lineC.split()
    '''



# Make a root file 
f = TFile('StatesCompare.root','RECREATE')
gStyle.SetOptStat(111111) #works on Draw() only 

h_Good_C_hits  = TH1F("h_Good_C_hits", "Good_C_hits [cm]", 100, -0.5, 0.5)
h_Bad_C_hits  = TH1F("h_Bad_C_hits", "Bad_C_hits [cm]", 100, -0.5, 0.5)
h_Good_Bad_C_hits  = TH1F("h_Good_Bad_C_hits ", "$\Delta$ Good-Bad_C_hits [cm]", 100, -0.000004, 0.000004)

print len(Good_C_hits)
BOOOcounter = 0
for n in range(0, len(Good_C_hits)): 
    h_Good_C_hits.Fill(float(Good_C_hits[n]))
    h_Bad_C_hits.Fill(float(Bad_C_hits[n]))    
    h_Good_Bad_C_hits.Fill(float(Good_C_hits[n]) - float(Bad_C_hits[n]))
    if (float(Good_C_hits[n]) - float(Bad_C_hits[n]) != 0):
        print "Hits are not the same BOOOOOOOO! ", Good_C_hits[n], " vs ", Bad_C_hits[n], " at ", n
        BOOOcounter=BOOOcounter+1

print "Total BOOOs are at ", BOOOcounter


'''
h_y0.Draw()
print "Hit Enter to continue"
test = raw_input()
'''

f.Write()
f.Close()








