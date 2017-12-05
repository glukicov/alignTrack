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
import ROOT
ROOT.gROOT.Macro('rootlogon.C')
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
        derLC_2.append(float(number[2]))
        ngl.append(number[3])
        derGL_1.append(float(number[4]))
        label_1.append(number[5])
        hits.append(float(number[6]))
        errors.append(float(number[7]))
           
  
# Make a root file 
f = TFile('MP2Tracker.root','RECREATE')

h_nlc  = TH1I("h_nlc", "NLC", 500, 0, 4)
h_ngl  = TH1I("h_ngl", "NGL", 500, 0, 2)
h_label_1  = TH1I("h_label_1", "label_1", 100, 0, 8)
h_derLc1_P1  = TH1F("h_derLc1_P1", "lc1 zoom in at +1", 1000, 0.9986, 1.0002)
h_derLc1_M1  = TH1F("h_derLc1_M1", "lc1 zoom in at -1", 1000, -1.0001, -0.9998)
h_derLc1  = TH1F("h_derLc1", "lc1", 500, -1.4, 1.4)
h_derLc2  = TH1F("h_derLc2", "lc2", 500, -65, 65)
h_derLc2_Zoom  = TH1F("h_derLc2_Zoom", "lc2 Zoom ", 150, 41.098, 41.106)
h_derGl1  = TH1F("h_derGl1", "gl1", 500, -1.5, 1.5)
h_hits  = TH1F("h_hits", "hits [cm]", 279, -0.08, 0.08)
h_errors  = TH1F("h_erros", "errors [cm]", 100, 0.01, 0.02)

for n in range(0, len(nlc)): 
    h_nlc.Fill(int(nlc[n]))
    h_ngl.Fill(int(ngl[n]))
    h_label_1.Fill(int(label_1[n]))
    h_derLc1.Fill(float(derLC_1[n]))
    h_derLc1_M1.Fill(float(derLC_1[n]))
    h_derLc1_P1.Fill(float(derLC_1[n]))
    h_derLc2.Fill(float(derLC_2[n]))
    h_derLc2_Zoom.Fill(float(derLC_2[n]))
    h_derGl1.Fill(float(derGL_1[n]))
    h_hits.Fill(float(hits[n]))
    h_erros.Fill(float(errors[n]))


gStyle.SetOptStat(111111)
c1 = TCanvas("c1","c1",1700,1500)
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

c2 = TCanvas("c2","c2",1700,1500)
c2.Divide(2,2)
c2.cd(1) 
h_label_1.Draw()
c2.cd(2) 
h_hits.Draw()
c2.cd(3) 
h_erros.Draw()
c2.cd(4)
h_derLc2.Draw()
c2.Print("c2_mp2.png")

f.Write()
f.Close()

print "ROOT file ready: root MP2Tracker.root"
print "Canvases Printed to files: c1_mp2.png and c2_mp2.png"






