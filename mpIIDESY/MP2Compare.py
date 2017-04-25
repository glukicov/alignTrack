#!/usr/bin/env python

import os
import string
import time
import decimal
from ROOT import *
#import numpy as np # TODO 

'''
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
    gStyle.SetOptFit(1111)
  
    gStyle.SetMarkerStyle(21)
    gStyle.SetMarkerSize(0.8)
    gStyle.SetMarkerColor(ROOT.kBlue)
    '''

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
C_projX=[]
C_hitsX=[]
C_x=[]
C_SdevX=[]

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
F_projX=[]
F_hitsX=[]
F_x=[]
F_SdevX=[]


#Read file and store in lists 
with open("../C_mpIIDESY/C_Mp2debug_mp2.txt") as f:
    #next(f)
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
        C_projX.append(number[12])
        C_hitsX.append(number[13])
           
        
   
with open("../F_mpIIDESY/F_mp2test2_mp2_debug.txt") as f:
    #next(f)
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
        F_projX.append(number[12])
        F_hitsX.append(number[13])

C_xl=[]
C_xs=[]
C_projX2=[]
C_yl=[]
C_ys=[]
C_projY=[]
C_gaus=[]
C_res=[]

F_xl=[]
F_xs=[]
F_projX2=[]
F_yl=[]
F_ys=[]
F_projY=[]
F_gaus=[]
F_res=[]

with open("../C_mpIIDESY/C_Mp2debug_off.txt") as f:
    #next(f)
    for lineC in f:  
        number = lineC.split()
        C_xl.append(number[0])
        C_xs.append(number[1])
        C_projX2.append(number[2])
        C_yl.append(number[3])
        C_ys.append(number[4])
        C_projY.append(number[5])
        C_gaus.append(number[6])
        C_res.append(number[7])
        C_x.append(number[8])
        C_SdevX.append(number[9])

        
with open("../F_mpIIDESY/F_mp2test2_off_debug.txt") as f:
    #next(f)
    for lineF in f:  
        number = lineF.split()
        F_xl.append(number[0])
        F_xs.append(number[1])
        F_projX2.append(number[2])
        F_yl.append(number[3])
        F_ys.append(number[4])
        F_projY.append(number[5])
        F_gaus.append(number[6])
        F_res.append(number[7])
        F_x.append(number[8])
        F_SdevX.append(number[9])



# Sanity Checks for labels
list_lenght=len(C_hits)
list_lenght2=len(C_xl) 
#Root_StyleSetup()
# Make a root file 
f = TFile('MP2Compare.root','RECREATE')
h_nlc  = TH1I("h_nlc", "$\Delta$ NLC", 100, -0.0000001, 0.0000001)
h_ngl  = TH1I("h_ngl", "$\Delta$ NGL", 100, -0.0000001, 0.0000001)
h_label_1  = TH1I("h_label_1", "$\Delta$ label 1", 100, -0.0000001, 0.0000001)
h_label_2  = TH1I("h_label_2", "$\Delta$ label 1", 100, -0.0000001, 0.0000001)

h_lc1  = TH1F("h_lc1", "$\Delta$ lc1", 300, -0.0000001, 0.0000001)
h_lc2  = TH1F("h_lc2", "$\Delta$ lc2", 100, -0.0000001, 0.0000001)
h_lc3  = TH1F("h_lc3", "$\Delta$ lc3", 300, -0.00001, 0.0001)
h_lc4  = TH1F("h_lc4", "$\Delta$ lc4", 300, -0.00001, 0.000001)
h_gl1  = TH1F("h_gl1", "$\Delta$ gl1", 300, -0.0000001, 0.0000001)
h_gl2  = TH1F("h_gl2", "$\Delta$ gl1", 100, -0.0000001, 0.0000001)

h_hits  = TH1F("h_hits", "$\Delta$ hits [cm]", 100, -0.02, 0.005)
h_errors  = TH1F("h_erros", "$\Delta$ errors [cm]", 100, -0.0000001, 0.0000001)

h_projX =  TH1F("h_projX", "$\Delta$ projX [cm]", 300, -0.0000001, 0.0000001)
h_hitsX =  TH1F("h_hitsX", "$\Delta$ hitsX [cm]", 100, -0.0000001, 0.0000001)



h_xl =  TH1F("h_xl", "$\Delta$ xl [cm]", 300, -0.05, 0.05)
h_xs =  TH1F("h_xs", "$\Delta$ xs [cm]", 100, -0.000005, 0.000005,)
h_projX2 =  TH1F("h_projX2", "$\Delta$ projX2 [cm]", 300, -0.000003, 0.000003)

h_yl =  TH1F("h_yl", "$\Delta$ yl [cm]", 300, -0.05, 0.05)
h_ys =  TH1F("h_ys", "$\Delta$ ys [cm]", 100, -0.000003, 0.000003)
h_projY =  TH1F("h_projY", "$\Delta$ projY [cm]", 300, -0.0000001, 0.0000001)

h_gaus =  TH1F("h_gaus", "$\Delta$ gaus [cm]", 300, -0.0000001, 0.0000001)
h_res =  TH1F("h_res", "$\Delta$ res [cm]", 300, -0.0000001, 0.0000001)

h_x =  TH1F("h_x", "$\Delta$ x [cm]", 300, -0.0004, 0.0004)
h_SdevX =  TH1F("h_SdevX", "$\Delta$ SdevX [cm]", 300, -0.005, 0.02)


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

    h_projX.Fill(float(C_projX[n])-float(F_projX[n]))
    h_hitsX.Fill(float(C_hitsX[n])-float(F_hitsX[n]))


for n in range(0, list_lenght2-1):
    h_xl.Fill(float(C_xl[n])-float(F_xl[n]))
    h_xs.Fill(float(C_xs[n])-float(F_xs[n]))
    h_projX2.Fill(float(C_projX2[n])-float(F_projX2[n]))

    h_yl.Fill(float(C_yl[n])-float(F_yl[n]))
    h_ys.Fill(float(C_ys[n])-float(F_ys[n]))
    h_projY.Fill(float(C_projY[n])-float(F_projY[n]))

    h_gaus.Fill(float(C_gaus[n])-float(F_gaus[n]))
    h_res.Fill(float(C_res[n])-float(F_res[n]))

    h_x.Fill(float(C_x[n])-float(F_x[n]))
    h_SdevX.Fill(float(C_SdevX[n])-float(F_SdevX[n]))

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

#f.SetOptStat(111111)
f.Write()
f.Close()








