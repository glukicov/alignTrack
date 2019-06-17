import sys
from ROOT import TFile, TStyle, TCanvas, gStyle, TF1, gROOT
import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines 
import argparse, sys
from math import log10, floor
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
from array import array
import subprocess
def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)

from ROOT import TH1F, TH2F, TF1, TCanvas, TFile, gStyle, TPaveText, TLegend, TProfile, TGraphErrors
from decimal import *
round_to = 3
getcontext().prec = round_to


gStyle.SetOptStat(0) 
gStyle.SetOptFit(0)
gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.033)
# gStyle.ForceStyle()

#Define constant paths and labels 
path =  "MomentumSlices/vertices/station"
TfileName = "gm2tracker_MomSlices_ana.root"
stationName = ["12", "18"]
topDir= "/Users/gleb/software/alignTrack/mpIIDESY/BK_AlignCurvature_2D/"
states = ["run15922_align_noVolumes", "sim_a=-0.5e-6", "sim_truth", "sim_a=1e-7", "sim_a=0.2e-6", "sim_a=0.3e-6", "sim_a=0.4e-6", "sim_a=0.5e-6", "sim_a=1e-6"]
stateN=len(states)
names = ["run15922 aligned", "sim_a=-0.5e-6", "sim_truth", "sim_a=1e-7", "sim_a=0.2e-6", "sim_a=0.3e-6", "sim_a=0.4e-6", "sim_a=0.5e-6", "sim_a=1e-6"]

#Containers to store histograms in orders as the names 

colors = [1, 2, 3 ,4 ,5, 6 ,7 ,8 ,9] #purple, green 
styles = [3001, 3002]
plotName = ["h_radialPos_vs_mom_timeCut", "h_radialPos_vs_mom", "h_radialPos_vs_mom", "h_radialPos_vs_mom", "h_radialPos_vs_mom", "h_radialPos_vs_mom", "h_radialPos_vs_mom", "h_radialPos_vs_mom", "h_radialPos_vs_mom"]
plotYtitle = ["Radial Beam Position"]
plotXtitle = ["P"]
plotTitle = [" "]
plotMean = [" <x>"]
unitsY = ["mm"]
unitsX = ["MeV"]

meanArray=[] # for the final FoM shift-nominal 

#Make new canvas for plots 
c = TCanvas("h_radialPos_vs_mom_timeCut", "h_radialPos_vs_mom_timeCut", 2800, 3800)
c.Divide(1,2)
#Keep legend, histots and TFiles in scope 
legendArray=[]
histArray=[]
fileArray=[]
profileArray=[]
graphArray=[]

i_total=0 # canvas id counter 
for i_station in range(0, len(stationName)):
    #canvas pad and legend per station 
    c.cd(i_total+1)
    legend =  TLegend(0.30,0.12,0.45,0.45)
    legendArray.append(legend) # stroe all to keep in scope 
    for i_plot in range(0, stateN):
    
        #Open TFiles
        fullPath = topDir+states[i_plot]+"/"+TfileName
        #print(fullPath)
        scrFile = TFile.Open(fullPath)
        fileArray.append(scrFile)

        #Get the TH2F 
       # print(path+stationName[i_station]+"/"+plotName[i_plot])
        plot = scrFile.Get(str(path+stationName[i_station]+"/"+plotName[i_plot]))
        histArray.append(plot)

        #Get the normalisation scale
        correction = plot.GetMean(2)
        #print(correction)

        #Rebin2D(X,Y)
        # plot.Rebin2D(1,1)

        #Get position stats along the y-axis 
        mean=plot.GetMean(2)
        mean_error = plot.GetMeanError(2)
        sd=plot.GetRMS(2)
        sd_error=plot.GetRMSError(2)
        meanArray.append(mean)

        #plot tools 
        binN_Y=plot.GetYaxis().GetBinWidth(1)
        binN_X=plot.GetXaxis().GetBinWidth(1)
        
        #make a profile 
        profile = plot.ProfileX("profile_"+str(i_plot)+str(i_station), 1, -1, "")
        profileArray.append(profile)
        profile.SetMarkerStyle(20)

        # create a new TGrpah for data with the same number of bins as the profile
        # only for non-zero bins
        bin_number_X = profile.GetNbinsX()
        x=array( 'f', [])
        y=array( 'f', [])
        ex=array( 'f', [])
        ey=array( 'f', [])

        # for i_bin_x in range(hBinMin, hBinMax+1):
        binNumber_actual = 0
        for i_bin_x in range(0, bin_number_X):
            bin_content = profile.GetBinContent(i_bin_x)
            if (bin_content == 0):
                continue
            new_bin_content = bin_content - correction
            y.append(new_bin_content)
            x.append(profile.GetBinCenter(i_bin_x))
            ex.append(profile.GetBinWidth(i_bin_x))
            ey.append(profile.GetBinError(i_bin_x))
            binNumber_actual+=1

        gr = TGraphErrors(binNumber_actual, x, y ,ex, ey)
        gr.SetName("TGraphErrors_"+str(i_plot)+str(i_station))
        graphArray.append(gr)

        #new profile y-axis tools 
        profile.GetYaxis().SetTitle(plotYtitle[0]+"/ " +str(binN_Y)+" "+ unitsY[0])            
        profile.GetYaxis().SetTitleSize(0.045)            
        profile.GetYaxis().SetLabelSize(0.040)            
        profile.GetXaxis().SetTitleSize(0.045)            
        profile.GetXaxis().SetLabelSize(0.040)            
        profile.GetYaxis().CenterTitle()
        profile.GetYaxis().SetTitleOffset(0.45)
        profile.GetYaxis().SetRangeUser(-25.0, 11.0) 
        profile.SetTitle("S"+stationName[i_station])
        profile.GetXaxis().SetTitle(plotXtitle[0]+"/ "+str(round(binN_X,3))+" "+ unitsX[0])
        profile.GetXaxis().CenterTitle()

        #new profile y-axis tools (in case tgrpah is drawn first)
        gr.GetYaxis().SetTitle(plotYtitle[0]+"/ " +str(binN_Y)+" "+ unitsY[0])            
        gr.GetYaxis().SetTitleSize(0.045)            
        gr.GetYaxis().SetLabelSize(0.040)            
        gr.GetXaxis().SetTitleSize(0.045)            
        gr.GetXaxis().SetLabelSize(0.040)            
        gr.GetYaxis().CenterTitle()
        gr.GetYaxis().SetTitleOffset(0.45)
        gr.GetYaxis().SetRangeUser(-25.0, 11.0) 
        gr.SetTitle("S"+stationName[i_station])
        gr.GetXaxis().SetTitle(plotXtitle[0]+"/ "+str(round(binN_X,3))+" "+ unitsX[0])
        gr.GetXaxis().CenterTitle()
        
        gr.SetMarkerColor( colors[i_plot] )
        gr.SetLineColor( colors[i_plot] )
        #profile.SetMarkerColor(colors[i_plot])  # profile 
        #profile.SetLineColor(colors[i_plot])  # profile 
        if (i_plot == 0):  
            #profile.Draw("E")  # profile  
            gr.Draw("AP")
           
        else:
            # profile.Draw("E same")  # profile  
            gr.Draw("AP same")

        #fill legend once per state      
        # legenValue1 = str(names[i_plot])+plotMean[0]+": "+str(round(mean, round_to))+" #pm "+str(round(mean_error,round_to)) # profile 
        legenValue1 = str(names[i_plot])
        legend.AddEntry(gr, str(legenValue1))
        
    #draw legend once per canvas 
    legend.Draw("same")
    i_total+=1
    meanArray=[]
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)

c.Draw()
c.Print("h_radialPos_vs_mom_timeCut.png")
