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

parser = argparse.ArgumentParser(description='arguments')
parser.add_argument('-m', '--mode', type=str)
parser.add_argument('-c', '--coord', type=str, default="radial")
parser.add_argument('-xmin', '--x_min', type=float, default=300)
parser.add_argument('-xmax', '--x_max', type=float, default=3000)
parser.add_argument('-ymin', '--y_min', type=float, default=-25)
parser.add_argument('-ymax', '--y_max', type=float, default=11)
args = parser.parse_args()

gROOT.SetBatch(1)

mode = args.mode
if not (mode == "profile" or mode== "graph"):
    print("Specify 'profile' or 'graph' as --mode=")
    sys.exit()

coord = args.coord
if not (coord == "radial" or coord== "vertical"):
    print("Specify 'radial' or 'vertical' as --coord=")
    sys.exit()

y_min = args.y_min
y_max = args.y_max
x_min = args.x_min
x_max = args.x_max


gStyle.SetOptStat(0) 
gStyle.SetOptFit(0)
gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.033)
# gStyle.ForceStyle()

#Define constant paths and labels 
path =  "MomentumSlices/vertices/station"
TfileName = "gm2tracker_MomSlices_ana.root"
stationName = ["12", "18"]
# topDir= "/Users/gleb/software/alignTrack/mpIIDESY/AlignCurvature_2D/" 
states = ["run15922-15923_align", "sim_a=-0.5e-6", "sim_truth", "sim_a=1e-7", "sim_a=0.2e-6", "sim_a=0.3e-6", "sim_a=0.4e-6", "sim_a=0.5e-6", "sim_a=1e-6"]
stateN=len(states)
names = ["runs15922-15923 aligned", "sim_a=-0.5e-6", "sim_a=truth", "sim_a=0.1e-6", "sim_a=0.2e-6", "sim_a=0.3e-6", "sim_a=0.4e-6", "sim_a=0.5e-6", "sim_a=1.0e-6"]

#Containers to store histograms in orders as the names 

colors = [1, 2, 3 ,4 ,5, 6 ,7 ,41 ,9] #purple, green 
styles = [3001, 3002]
plotName = ["h_"+coord+"Pos_vs_mom_timeCut", "h_"+coord+"Pos_vs_mom"]
plotYtitle = [coord+" beam position"]
plotXtitle = ["P"]
plotTitle = [" "]
plotMean = ["<"+coord+">"]
unitsY = ["mm"]
unitsX = ["MeV"]


#Make new canvas for plots 
w = 1400
h = 600
c = TCanvas("h_radialPos_vs_mom_timeCut", "h_radialPos_vs_mom_timeCut", w, h)
c.SetWindowSize(w + (w - c.GetWw()), h + (h - c.GetWh()))
c.Divide(2,1)
#Keep legend, histots and TFiles in scope 
legendArray=[]
histArray=[]
fileArray=[]
profileArray=[]
graphArray=[]
meanArray=[] # for the final FoM shift-nominal 

i_total=0 # canvas id counter 
for i_station in range(0, len(stationName)):
    #canvas pad and legend per station 
    c.cd(i_total+1)
    legend =  TLegend(0.30,0.12,0.45,0.45)
    legendArray.append(legend) # stroe all to keep in scope 
    for i_plot in range(0, stateN):
        
        tmp_holder = [] # takes in profile or graph for the inner loop only 
        
        #Open TFiles
        fullPath = states[i_plot]+"/"+TfileName
        #print(fullPath)
        scrFile = TFile.Open(fullPath)
        fileArray.append(scrFile)

        #Get the TH2F 
       # print(path+stationName[i_station]+"/"+plotName[i_plot])
        if (i_plot == 0):
            plotNameAfterCut = plotName[0]
        else:
            plotNameAfterCut = plotName[1]
        plot = scrFile.Get(str(path+stationName[i_station]+"/"+plotNameAfterCut))
        histArray.append(plot)

       # print("before: ", plot.GetNbinsX())
        # Rebin2D(X,Y)
        # plot.Rebin2D(4,1)
        #print("after: ", plot.GetNbinsX())

        #Get the normalisation scale
        correction = plot.GetMean(2) 
        # print(correction)

        #Get position stats along the y-axis 
        mean=plot.GetMean(2)
        mean_error = plot.GetMeanError(2)
        sd=plot.GetRMS(2)
        sd_error=plot.GetRMSError(2)
        meanArray.append(mean)

        #plot tools 
        binN_Y=plot.GetYaxis().GetBinWidth(1)
        binN_X=plot.GetXaxis().GetBinWidth(1)
        
        #make a profileX  
        profile = plot.ProfileX("profile_"+str(i_plot)+str(i_station), 1, -1, "")
        profileArray.append(profile)
        
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
            ex.append(profile.GetBinWidth(i_bin_x)*0.5)
            ey.append(profile.GetBinError(i_bin_x))
            binNumber_actual+=1

        gr = TGraphErrors(binNumber_actual, x, y ,ex, ey)
        gr.SetName("TGraphErrors_"+str(i_plot)+str(i_station))
        graphArray.append(gr)

        #plot depending on the mode 
        tmp_holder.append(profile)
        tmp_holder.append(gr)
        holder_i = None
        if (mode == "profile"):
            holder_i = 0 
        if (mode == "graph"):
            holder_i = 1 
                
        #new profile/graph y-axis tools (in case tgrpah is drawn first)
        tmp_holder[holder_i].GetYaxis().SetTitle(plotYtitle[0]+"/ " +str(binN_Y)+" "+ unitsY[0])            
        tmp_holder[holder_i].GetYaxis().SetTitleSize(0.045)            
        tmp_holder[holder_i].GetYaxis().SetLabelSize(0.040)            
        tmp_holder[holder_i].GetXaxis().SetTitleSize(0.045)            
        tmp_holder[holder_i].GetXaxis().SetLabelSize(0.040)            
        tmp_holder[holder_i].GetYaxis().CenterTitle()
        tmp_holder[holder_i].GetYaxis().SetTitleOffset(0.8)
        tmp_holder[holder_i].GetYaxis().SetRangeUser(y_min, y_max) 
        tmp_holder[holder_i].GetXaxis().SetRangeUser(x_min, x_max) 
        tmp_holder[holder_i].SetTitle("S"+stationName[i_station])
        tmp_holder[holder_i].GetXaxis().SetTitle(plotXtitle[0]+"/ "+str(round(binN_X,3))+" "+ unitsX[0])
        tmp_holder[holder_i].GetXaxis().CenterTitle()
        tmp_holder[holder_i].SetMarkerColor( colors[i_plot] )
        tmp_holder[holder_i].SetLineColor( colors[i_plot] )
        tmp_holder[holder_i].SetMarkerStyle(20)

        if (i_plot == 0):  
            if (mode == "profile"):
                profile.Draw("E")  # profile
            if (mode == "graph"):
                gr.Draw("AP")
                #gr.SetMarkerSize(2)
                gr.SetMarkerStyle(20)
        else:
            if (mode == "profile"):
                profile.Draw("E same")  # profile  
            if (mode == "graph"):
                gr.Draw("P same")
            

        #fill legend once per state      
        if (mode == "profile"):
            legenValue1 = str(names[i_plot])+" "+plotMean[0]+": "+str(round(mean, round_to))+" #pm "+str(round(mean_error,round_to)) # profile 
        if (mode == "graph"):
            legenValue1 = str(names[i_plot])
        legend.AddEntry(tmp_holder[holder_i], str(legenValue1))
        
    #draw legend once per canvas 
    legend.Draw("same")
    i_total+=1
    meanArray=[]
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)

c.Draw()
c.Print(coord+"_"+mode+"_"+"Pos_vs_mom_timeCut.png")
c.SaveAs(coord+"_"+mode+"_"+"Pos_vs_mom_timeCut.C")