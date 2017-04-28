#!/usr/bin/env python

import os
import string
import time
import decimal
from ROOT import *
#import numpy as np # TODO 
    
C_labels=[]
F_labels=[]
C_results=[]
F_results=[]
C_errors=[]
F_errors=[]

#Read file and store in lists 
with open("../C_mpIIDESY/millepede.res") as f:
    next(f)
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number_str = line.split()    
        C_labels.append(number_str[0])
        C_results.append(number_str[1])
        C_errors.append(number_str[len(number_str)-1])
   
with open("millepede.res") as f:
    next(f)
    for line in f:  
        number_str = line.split()
        F_labels.append(number_str[0])
        F_results.append(number_str[1])
        F_errors.append(number_str[len(number_str)-1])


# Sanity Checks for labels
n = 0  # list element counter 
errorCounterLabels = 0
errorCounterResults = 0
errorCounterErrors = 0
CF_diff_results=[]
CF_diff_errors=[]
zeroCounter = 0

for item in C_labels:
    if (float(item) - float(F_labels[n]) != 0):
        errorCounterLabels = errorCounterLabels + 1
    if (float(C_results[n]) - float(F_results[n]) != 0):
        errorCounterResults = errorCounterResults + 1
    if (float(C_errors[n]) - float(F_errors[n]) != 0):
        errorCounterErrors = errorCounterErrors + 1

    if (float(C_results[n])==0.0 and float(F_results[n])==0.0):
        zeroCounter=zeroCounter+1 
    
    CF_diff_results.append(float(C_results[n])-float(F_results[n]))
    CF_diff_errors.append(float(C_errors[n])-float(F_errors[n]))

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
h_res  = TH1F("h_res", "C-F alignment results [cm]", 100, -0.0000002, 0.0000002)
 
for item in CF_diff_results: 
    h_res.Fill(item)

h_res.Draw()

print "Hit Enter to continue"
test = raw_input()


#TODO gROOT set opst stats etc. 

f.Write()
f.Close()








