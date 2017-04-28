#!/usr/bin/env python

import os
import string
import time
import decimal
from ROOT import *

Good_C_hits=[]
Bad_C_hits=[]

Good_C_rand4=[]
Bad_C_rand4=[]

Good_F_hits=[]
Bad_F_hits=[]

Good_F_rand4=[]
Bad_F_rand4=[]


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

h_Good_C_hits  = TH1F("h_Good_C_hits", "Good_C_hits [cm]", 300, -0.1, 0.1)
h_Bad_C_hits  = TH1F("h_Bad_C_hits", "Bad_C_hits [cm]", 300, -0.1, 0.1)
h_Good_Bad_C_hits  = TH1F("h_Good_Bad_C_hits ", "$\Delta$ Good-Bad_C_hits [cm]", 300, -0.000001, 0.000001)

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








