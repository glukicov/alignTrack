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

C_projX2=[]
C_projY=[]
C_gaus=[]
C_res=[]


F_projX2=[]
F_projY=[]
F_gaus=[]
F_res=[]

with open("../C_mpIIDESY/C_Mp2debug_off.txt") as f:
    #next(f)
    for lineC in f:  
        number = lineC.split()
        C_projX2.append(number[0])
        C_projY.append(number[1])
        C_gaus.append(number[2])
        C_res.append(number[3])
        C_SdevX.append(number[4])

        
with open("../F_mpIIDESY/F_mp2test2_off_debug.txt") as f:
    #next(f)
    for lineF in f:  
        number = lineF.split()
        F_projX2.append(number[0])
        F_projY.append(number[1])
        F_gaus.append(number[2])
        F_res.append(number[3])
        F_SdevX.append(number[4])


C_x0=[]
C_y0=[]
C_x1=[]
C_y1=[]
C_xslope=[]
C_yslope=[]
C_xl=[]
C_yl=[]
C_xs=[]
C_ys=[]
C_x=[]
C_y=[]
C_rand1=[]
C_rand2=[]
C_rand3=[]
C_rand4=[]
C_layerSize=[]
C_distance=[]
C_dx=[]
C_dy=[]

F_x0=[]
F_y0=[]
F_x1=[]
F_y1=[]
F_xslope=[]
F_yslope=[]
F_xl=[]
F_yl=[]
F_xs=[]
F_ys=[]
F_x=[]
F_y=[]
F_rand1=[]
F_rand2=[]
F_rand3=[]
F_rand4=[]
F_layerSize=[]
F_distance=[]
F_dy=[]
F_dx=[]

with open("../C_mpIIDESY/C_Mp2debug_calc.txt") as f:
    #next(f)
    for lineC in f:
        number = lineC.split()
        C_x0.append(number[0])
        C_y0.append(number[1])
        C_x1.append(number[2])
        C_y1.append(number[3])
        C_xslope.append(number[4])
        C_yslope.append(number[5])
        C_rand1.append(number[6])
        C_rand2.append(number[7])
        C_rand3.append(number[8])
        C_rand4.append(number[9])
        C_layerSize.append(number[10])
        C_distance.append(number[11])
        

with open("../F_mpIIDESY/F_mp2test2_calc_debug.txt") as f:
    #next(f)
    for lineF in f:  
        number = lineF.split()
        F_x0.append(number[0])
        F_y0.append(number[1])
        F_x1.append(number[2])
        F_y1.append(number[3])
        F_xslope.append(number[4])
        F_yslope.append(number[5])
        F_rand1.append(number[6])
        F_rand2.append(number[7])
        F_rand3.append(number[8])
        F_rand4.append(number[9])
        F_layerSize.append(number[10])
        F_distance.append(number[11])
    

with open("../C_mpIIDESY/C_Mp2debug_MC.txt") as f:
    #next(f)
    for lineC in f:
        number = lineC.split()
        C_xl.append(number[0])
        C_yl.append(number[1])
        C_xs.append(number[2])
        C_ys.append(number[3])
        C_dx.append(number[4])
        C_dy.append(number[5])
        C_x.append(number[6])
        C_y.append(number[7])

with open("../F_mpIIDESY/F_mp2test2_MC_debug.txt") as f:
    #next(f)
    for lineC in f:
        number = lineC.split()
        F_xl.append(number[0])
        F_yl.append(number[1])
        F_xs.append(number[2])
        F_ys.append(number[3])
        F_dx.append(number[4])
        F_dy.append(number[5])
        F_x.append(number[6])
        F_y.append(number[7])


# Sanity Checks for labels
list_lenght=len(C_nlc)
list_lenght2=len(C_res) 
list_lenght3=len(C_rand1)

list_lenght_C=len(C_rand1)
list_lenght_F=len(F_rand1)

list_lenghtdx=len(C_dx)

if (list_lenght_C != list_lenght_F):
    print "Rand lists are not equal!"

# Make a root file 
f = TFile('MP2Compare.root','RECREATE')
gStyle.SetOptStat(111111) #works on Draw() only 

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

h_hits  = TH1F("h_hits", "$\Delta$ hits [cm]", 100, -0.000005, 0.000005)
h_errors  = TH1F("h_erros", "$\Delta$ errors [cm]", 100, -0.0000001, 0.0000001)

h_projX =  TH1F("h_projX", "$\Delta$ projX [cm]", 300, -0.0000001, 0.0000001)
h_hitsX =  TH1F("h_hitsX", "$\Delta$ hitsX [cm]", 100, -0.0000001, 0.0000001)


h_projX2 =  TH1F("h_projX2", "$\Delta$ projX2 [cm]", 300, -0.000003, 0.000003)
h_projY =  TH1F("h_projY", "$\Delta$ projY [cm]", 300, -0.0000001, 0.0000001)

h_gaus =  TH1F("h_gaus", "$\Delta$ gaus [cm]", 300, -0.0000001, 0.0000001)
h_res =  TH1F("h_res", "$\Delta$ res [cm]", 300, -0.0000001, 0.0000001)

h_SdevX =  TH1F("h_SdevX", "$\Delta$ SdevX [cm]", 300, -0.005, 0.02)


h_xl =  TH1F("h_xl", "$\Delta$ xl [cm]", 50, -0.000015, 0.000015)
h_xs =  TH1F("h_xs", "$\Delta$ xs [cm]", 50, -0.000015, 0.000015)
h_x =  TH1F("h_x", "$\Delta$ x [cm]", 50, -0.00002, 0.00002)
h_y =  TH1F("h_y", "$\Delta$ y [cm]", 50, -0.00002, 0.00002)
h_yl =  TH1F("h_yl", "$\Delta$ yl [cm]", 50, -0.00002, 0.00002)
h_ys =  TH1F("h_ys", "$\Delta$ ys [cm]", 50, -0.000015, 0.000015)
h_dx =  TH1F("h_dx", "$\Delta$ yl [cm]", 50, -0.00002, 0.00002)
h_dy =  TH1F("h_dy", "$\Delta$ ys [cm]", 50, -0.000015, 0.000015)

h_xslope =  TH1F("h_xslope", "$\Delta$ xslope [cm]", 50, -0.000002, 0.000002)
h_yslope =  TH1F("h_yslope", "$\Delta$ yslope [cm]", 50, -0.000002, 0.000002)
h_x1 =  TH1F("h_x1", "$\Delta$ x1 [cm]", 50, -0.000002, 0.000002)
h_y1 =  TH1F("h_y1", "$\Delta$ y1 [cm]", 50, -0.000002, 0.000002)
h_x0 =  TH1F("h_x0", "$\Delta$ x0 [cm]", 50, -0.000002, 0.000002)
h_y0 =  TH1F("h_y0", "$\Delta$ y0 [cm]", 50, -0.000002, 0.000002)

h_rand1 =  TH1F("h_rand1", "$\Delta$ rand1 [cm]", 50, -0.000002, 0.000002)
h_rand2 =  TH1F("h_rand2", "$\Delta$ rand2 [cm]", 50, -0.000002, 0.000002)
h_rand3 =  TH1F("h_rand3", "$\Delta$ rand3 [cm]", 50, -0.000002, 0.000002)
h_rand4 =  TH1F("h_rand4", "$\Delta$ rand4 [cm]", 50, -0.000002, 0.000002)
h_layerSize =  TH1F("h_layerSize", "$\Delta$ layerSize [cm]", 50, -0.000001, 0.000001)
h_distance =  TH1F("h_distance", "$\Delta$ layerSize [cm]", 50, -0.000002, 0.000002)

print len(C_nlc)
for n in range(0, list_lenght): 
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


for n in range(0, list_lenght2):
    h_projX2.Fill(float(C_projX2[n])-float(F_projX2[n]))

    h_projY.Fill(float(C_projY[n])-float(F_projY[n]))

    h_gaus.Fill(float(C_gaus[n])-float(F_gaus[n]))
    h_res.Fill(float(C_res[n])-float(F_res[n]))

    h_SdevX.Fill(float(C_SdevX[n])-float(F_SdevX[n]))

for n in range(0, list_lenght_C):


    h_xslope.Fill(float(C_xslope[n])-float(F_xslope[n]))
    h_yslope.Fill(float(C_yslope[n])-float(F_yslope[n]))
    h_y1.Fill(float(C_y1[n])-float(F_y1[n]))
    h_x1.Fill(float(C_x1[n])-float(F_x1[n]))
    h_x0.Fill(float(C_x0[n])-float(F_x0[n]))
    h_y0.Fill(float(C_y0[n])-float(F_y0[n]))

    if ((float(C_x0[n])-float(F_x0[n])) != 0):
        print "x0 not equal", C_x0[n], F_x0[n], " at ", n, "from ", C_rand1[n], " and ", F_rand1[n]

    h_rand1.Fill(float(C_rand1[n])-float(F_rand1[n]))
    h_rand2.Fill(float(C_rand2[n])-float(F_rand2[n]))
    h_rand3.Fill(float(C_rand3[n])-float(F_rand3[n]))
    h_rand4.Fill(float(C_rand4[n])-float(F_rand4[n]))
    h_layerSize.Fill(float(C_layerSize[n])-float(F_layerSize[n]))
    h_distance.Fill(float(C_distance[n])-float(F_distance[n]))
    #print n

for n in range(0, list_lenghtdx):
    h_xl.Fill(float(C_xl[n])-float(F_xl[n]))
    h_yl.Fill(float(C_yl[n])-float(F_yl[n]))
    h_xs.Fill(float(C_xs[n])-float(F_xs[n]))
    h_ys.Fill(float(C_ys[n])-float(F_ys[n]))
    h_dx.Fill(float(C_dx[n])-float(F_dx[n]))
    h_dy.Fill(float(C_dy[n])-float(F_dy[n]))
    h_x.Fill(float(C_x[n])-float(F_x[n]))
    h_y.Fill(float(C_y[n])-float(F_y[n]))

'''
    if ( (float(C_rand1[n])-float(F_rand1[n]) != 0)):
        print "rand1 not equal", C_rand1[n], F_rand1[n], " at ", n
    if ( (float(C_rand2[n])-float(F_rand2[n]) != 0)):
        print "rand2 not equal", C_rand2[n], F_rand2[n], " at ", n
    if ( (float(C_rand3[n])-float(F_rand3[n]) != 0)):
        print "rand3 not equal", C_rand3[n], F_rand3[n], " at ", n
    if ( (float(C_rand4[n])-float(F_rand4[n]) != 0)):
        print "rand4 not equal", C_rand4[n], F_rand4[n], " at ", n
'''





'''
h_y0.Draw()
print "Hit Enter to continue"
test = raw_input()
'''

f.Write()
f.Close()








