#!/usr/bin/env python

import os
import string
import time
import decimal
from ROOT import *
#import numpy as np # TODO 


nlc=[]
ngl=[]
label_1=[]
label_2=[]
derLC_1 = []
derLC_2 = []
derLC_3 = []
derLC_4 = []
derGL_1 = []
derGL_2 = []
hits=[]
errors=[]
projX=[]
hitsX=[]
x=[]
SdevX=[]

#Read file and store in lists 
with open("Tracker_d_mille.txt") as f:
    #next(f)
    for lineC in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number = lineC.split()
        nlc.append(number[0])
        derLC_1.append(number[1])
        derLC_2.append(number[2])
        ngl.append(number[3])
        derGL_1.append(number[4])
        label_1.append(number[5])
        hits.append(number[6])
        errors.append(number[7])
           
  
# Make a root file 
f = TFile('MP2Tracker.root','RECREATE')

h_nlc  = TH1I("h_nlc", "NLC", 500, 0, 3)
h_ngl  = TH1I("h_ngl", "NGL", 500, 0, 2)
h_label_1  = TH1I("h_label_1", "$label 1", 100, 0, 10)
h_lc1  = TH1F("h_lc1", "lc1", 500, 0, 200)
h_lc2  = TH1F("h_lc2", "lc2", 1000, 0, 1200)
h_gl1  = TH1F("h_gl1", "gl1", 500, 0, 1000)
h_hits  = TH1F("h_hits", "hits [cm]", 100, 0, 1)
h_errors  = TH1F("h_erros", "errors [cm]", 100, 0, 0.1)

for n in range(0, len(nlc)): 
    h_nlc.Fill(int(nlc[n]))
    h_ngl.Fill(int(ngl[n]))
    h_label_1.Fill(int(label_1[n]))
    h_lc1.Fill(float(derLC_1[n]))
    h_lc2.Fill(float(derLC_2[n]))
    h_gl1.Fill(float(derGL_1[n]))
    h_hits.Fill(float(hits[n]))
    h_erros.Fill(float(errors[n]))


f.Write()
f.Close()

print "ROOT file ready: root MP2Tracker.root"






