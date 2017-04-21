#!/usr/bin/env python

import os
import string
import time
import decimal
from ROOT import *
#import numpy as np # TODO 
    
C_hits=[]
C_error=[]
C_derlc=[]
C_dergl=[]
C_nlc=[]
C_ngl=[]
C_errors=[]
C_labels=[]


F_hits=[]
F_error=[]
F_derlc=[]
F_dergl=[]
F_nlc=[]
F_ngl=[]
F_errors=[]
F_labels=[]

#Read file and store in lists 
with open("../C_mpIIDESY/C_Mp2debug_mp2.txt") as f:
    next(f)
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number = line.split()
        C_ngl.append(float(number[0]))
        C_nlc.append(float(number[1]))
        C_labels.append(float(number[2]))
        C_labels.append(float(number[3]))
        C_derlc.append(float(number[4]))
        C_derlc.append(float(number[5]))   
        C_derlc.append(float(number[6]))   
        C_derlc.append(float(number[7]))   
        C_dergl.append(float(number[8]))
        C_dergl.append(float(number[9]))
        C_hits.append(float(number[10]))
        C_errors.append(float(number[11]))   
        
   
with open("../F_mpIIDESY/F_mp2test2_mp2_debug.txt") as f:
    next(f)
    for line in f:  
        number_str = line.split()
        F_ngl.append(float(number[0]))
        F_nlc.append(float(number[1]))
        F_labels.append(float(number[2]))
        F_labels.append(float(number[3]))
        F_derlc.append(float(number[4]))
        F_derlc.append(float(number[5]))   
        F_derlc.append(float(number[6]))   
        F_derlc.append(float(number[7]))   
        F_dergl.append(float(number[8]))
        F_dergl.append(float(number[9]))
        F_hits.append(float(number[10]))
        F_errors.append(float(number[11])) 
        


# Sanity Checks for labels
n = 0  # list element counter 
nlc_counter = 0
ngl_counter = 0
derlc_counter = 0
CF_diff_results=[]
CF_diff_errors=[]

for item in C_labels:
    

    n=n+1  

if (errorCounterLabels == 0): 
    print "All labels agree!"
else:
    print "Total number of errors in labels is ", errorCounterLabels

if (errorCounterResults == 0): 
    print "All results agree!"
else:
    print "Total number of disagreements in labels is ", errorCounterResults, "out of ", len(CF_diff_results)
    print "There were ", zeroCounter, "'0' displacemnt results"
if (errorCounterErrors == 0): 
    print "All errors agree!"
else:
    print "Total number of errors in labels is ", errorCounterErrors 

# Make a root file 
f = TFile('TProfile.root','RECREATE')
h_res  = TH1F("h_res", "C-F alignment results [cm]", 100, -0.025, 0.017)
 
for item in CF_diff_results: 
    h_res.Fill(item)

h_res.Draw()

print "Hit Enter to continue"
test = raw_input()


#TODO gROOT set opst stats etc. 

f.Write()
f.Close()








