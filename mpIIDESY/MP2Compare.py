#!/usr/bin/env python

import os
import string
import time
import decimal
from ROOT import *
#import numpy as np # TODO 

def Root_StyleSetup():
    gROOT.Reset()
  
    gROOT.SetStyle("Plain")
    gROOT.ForceStyle()
  
    gStyle.SetCanvasColor(0)
    gStyle.SetPadTopMargin(.1)
    gStyle.SetPadBottomMargin(.1)
  
    gStyle.SetPadLeftMargin(.15)
    gStyle.SetPadRightMargin(.15)
    gStyle.SetPadBorderMode(0)
  
    gStyle.SetTitleYOffset(1.2)
    gStyle.SetTitleXOffset(1.2)
  
    gStyle.SetOptStat(111111)
    #gStyle.SetOptFit(1111)
  
    gStyle.SetMarkerStyle(21)
    gStyle.SetMarkerSize(0.8)
    gStyle.SetMarkerColor(ROOT.kBlue)

C_nlc=[]
C_ngl=[]
C_label_1=[]
C_label_2=[]
C_derLC_1 = []
C_derLC_2 = []
C_derLC_3 = []
C_derLC_4 = []
C_derGL_1 = []
C_derGL_2 = []
C_hits=[]
C_errors=[]

F_nlc=[]
F_ngl=[]
F_label_1=[]
F_label_2=[]
F_derLC_1 = []
F_derLC_2 = []
F_derLC_3 = []
F_derLC_4 = []
F_derGL_1 = []
F_derGL_2 = []
F_hits=[]
F_errors=[]


#Read file and store in lists 
with open("../C_mpIIDESY/C_Mp2debug_mp2.txt") as f:
    next(f)
    for lineC in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number = lineC.split()
        C_nlc.append(number[0])
        C_ngl.append(number[1])
        C_label_1.append(number[2])
        C_label_2.append(number[3])
        C_derLC_1.append(number[4])
        C_derLC_2.append(number[5])
        C_derLC_3.append(number[6])
        C_derLC_4.append(number[7])
        C_derGL_1.append(number[8])
        C_derGL_2.append(number[9])
        C_hits.append(number[10])
        C_errors.append(number[11])
           
        
   
with open("../F_mpIIDESY/F_mp2test2_mp2_debug.txt") as f:
    next(f)
    for lineF in f:  
        number = lineF.split()
        F_nlc.append(number[0])
        F_ngl.append(number[1])
        F_label_1.append(number[2])
        F_label_2.append(number[3])
        F_derLC_1.append(number[4])
        F_derLC_2.append(number[5])
        F_derLC_3.append(number[6])
        F_derLC_4.append(number[7])
        F_derGL_1.append(number[8])
        F_derGL_2.append(number[9])
        F_hits.append(number[10])
        F_errors.append(number[11])
        
        

# Sanity Checks for labels
list_lenght=len(C_hits) 
Root_StyleSetup()
# Make a root file 
f = TFile('MP2Compare.root','RECREATE')
h_nlc  = TH1I("h_nlc", "$\Delta$ NLC", 100, -1, 1)
h_ngl  = TH1I("h_ngl", "$\Delta$ NGL", 100, -1, 1)
h_label_1  = TH1I("h_label_1", "$\Delta$ label 1", 100, -1, 1)
h_label_2  = TH1I("h_label_2", "$\Delta$ label 1", 100, -1, 1)

h_lc1  = TH1F("h_lc1", "$\Delta$ lc1", 300, -0.0002, 0.0002)
h_lc2  = TH1F("h_lc2", "$\Delta$ lc2", 100, -0.005, 0.005)
h_lc3  = TH1F("h_lc3", "$\Delta$ lc3", 300, -0.0002, 0.0002)
h_lc4  = TH1F("h_lc4", "$\Delta$ lc4", 300, -0.0002, 0.0002)
h_gl1  = TH1F("h_gl1", "$\Delta$ gl1", 300, -0.0002, 0.0002)
h_gl2  = TH1F("h_gl2", "$\Delta$ gl1", 100, -0.005, 0.005)

h_hits  = TH1F("h_hits", "$\Delta$ hits [cm]", 100, -0.01, 0.005)
h_errors  = TH1F("h_erros", "$\Delta$ errors [cm]", 100, -0.000005, 0.000005)
 
for n in range(0, list_lenght-1): 
    h_nlc.Fill(int(C_nlc[n])-int(F_nlc[n]))
    h_ngl.Fill(int(C_ngl[n])-int(F_ngl[n]))
    h_label_1.Fill(int(C_label_1[n])-int(F_label_1[n]))
    h_label_2.Fill(int(C_label_2[n])-int(F_label_2[n]))

    h_lc1.Fill(float(C_derLC_1[n])-float(F_derLC_1[n]))
    h_lc2.Fill(float(C_derLC_2[n])-float(F_derLC_2[n]))
    h_lc3.Fill(float(C_derLC_3[n])-float(F_derLC_3[n]))
    h_lc4.Fill(float(C_derLC_4[n])-float(F_derLC_4[n]))

    h_gl1.Fill(float(C_derGL_1[n])-float(F_derGL_1[n]))
    h_gl2.Fill(float(C_derGL_2[n])-float(F_derGL_2[n]))

    h_hits.Fill(float(C_hits[n])-float(F_hits[n]))
    h_erros.Fill(float(C_errors[n])-float(F_errors[n]))

'''
h_hits.Draw()
print "Hit Enter to continue"
test = raw_input()
h_erros.Draw()
print "Hit Enter to continue"
test = raw_input()
'''

#TODO gROOT set opst stats etc. 

#Root_StyleSetup()

f.Write()
f.Close()








