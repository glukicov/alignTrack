#!/usr/bin/env python


####################################################################
# PEDE inputs plots for AlginTracker
#
# 
#
# Created: 16 May 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 16 May 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################

import os
import string
import time
import decimal
from ROOT import *


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
        derLC_1.append(float(number[1]))
        ngl.append(number[2])
        derGL_1.append(float(number[3]))
        label_1.append(number[4])
        hits.append(float(number[5]))
        errors.append(float(number[6]))
           
  
# Make a root file 
f = TFile('MP2Tracker.root','RECREATE')

h_nlc  = TH1I("h_nlc", "NLC", 500, 0, 2)
h_ngl  = TH1I("h_ngl", "NGL", 500, 0, 2)
h_label_1  = TH1I("h_label_1", "label_1", 100, 0, 8)
h_derLc1  = TH1F("h_derLc1", "lc1", 500, 0, 2)
h_derGl1  = TH1F("h_derGl1", "gl1", 500, 0, 2)
h_hits  = TH1F("h_hits", "hits [cm]", 400, -0.2, 0.2)
h_errors  = TH1F("h_erros", "errors [cm]", 100, 0.01, 0.02)

for n in range(0, len(nlc)): 
    h_nlc.Fill(int(nlc[n]))
    h_ngl.Fill(int(ngl[n]))
    h_label_1.Fill(int(label_1[n]))
    h_derLc1.Fill(float(derLC_1[n]))
    h_derGl1.Fill(float(derGL_1[n]))
    h_hits.Fill(float(hits[n]))
    h_erros.Fill(float(errors[n]))


c1 = TCanvas("c1","c1",700,900)
c1.Divide(2,2)
c1.cd(1) 
h_nlc.Draw()
c1.cd(2) 
h_derLc1.Draw()
c1.cd(3) 
h_ngl.Draw()
c1.cd(4) 
h_derGl1.Draw()
c1.Print("c1_mp2.png")

c2 = TCanvas("c2","c2",700,900)
c2.Divide(2,2)
c2.cd(1) 
h_label_1.Draw()
c2.cd(2) 
h_hits.Draw()
c2.cd(3) 
h_erros.Draw()
c2.Print("c2_mp2.png")


f.Write()
f.Close()

print "ROOT file ready: root MP2Tracker.root"
print "Canvases Printed to files: c1_mp2.png and c2_mp2.png"






