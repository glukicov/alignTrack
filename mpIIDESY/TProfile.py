#!/usr/bin/env python

import os
import string
import time
import decimal
from ROOT import *
#import numpy as np # TODO 
    
C=[]  #stores all numbers 
F=[]

#Read file and store in lists 
with open("C_Mp2debug_off.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number_str = line.split()    
        C.append(number_str[0])
   
with open("F_mp2test2_off_debug.txt") as f:
    for line in f:  
        number_str = line.split()
        F.append(number_str[0])


# Make a root file 
f = TFile('TProfile.root','RECREATE')
hprof  = TProfile("hprof","Profile of Normalised F-C difference wrt to C result",14,-1,16, 's');  #name title, NbinsX, xMin, xMax 
hprof_abs  = TProfile("hprof_abs","Profile of Absolute Normalised F-C difference wrt to C result",14,-1,16, 's');  #name title, NbinsX, xMin, xMax

i=int(0)  #to loop vales over 14 hits from the same line [included are also rejected "hits" so its ALWAYS 14] 
n=0  #to loop over Fortran 
for number in C: 
    i=i+int(1)
    y=(float(number)-float(F[n]))/float(number)
    hprof.Fill(i,y)
    hprof_abs.Fill(i,abs(y))
    n=n+1
    if (i==14):
      i=0

hprof.Draw()

print "Hit Enter to continue"
test = raw_input()

hprof_abs.Draw()

print "Hit Enter to continue"
test = raw_input()


f.Write()
f.Close()







