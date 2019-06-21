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
parser.add_argument('-ymax', '--y_max', type=float, default=25)
args = parser.parse_args()

gROOT.SetBatch(1)

mode = args.mode
if not (mode == "profile" or mode== "graph"):
    print("Specify 'profile' or 'graph' as --mode=")
    sys.exit()

coord = args.coord
if not (coord == "radial" or coord== "vertical" or coord== "vertex"):
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
# states = ["run15922-15924_align", "sim_truth_default", "sim_truth_cbo"]
states = ["run15922-15924_align", "sim_a=-0.5e-6", "sim_truth", "sim_a=1e-7", "sim_a=0.2e-6", "sim_a=0.3e-6", "sim_a=0.4e-6", "sim_a=0.5e-6", "sim_a=1e-6"]
stateN=len(states)
# names = ["data (aligned)", "sim default", "sim cbo"]
names = ["data (aligned)", "sim a=-0.5e-6", "sim a=truth  ", "sim a=0.1e-6 ", "sim a=0.2e-6 ", "sim a=0.3e-6 ", "sim a=0.4e-6 ", "sim a=0.5e-6 ", "sim a=1.0e-6 "]

#Containers to store histograms in orders as the names 

colors = [1, 2, 3 ,4 ,5, 6 ,7 ,41 ,9] #purple, green 
styles = [3001, 3002]
plotName = ["h_"+coord+"Pos_vs_mom_timeCut", "h_"+coord+"Pos_vs_mom"]
if (coord == "vertex"):
    plotName = ["h_radialPos_vs_vertexExtrapDist_timeCut", "h_radialPos_vs_vertexExtrapDist"]
plotYtitle = [coord+" beam position"]
plotXtitle = ["P"]
unitsX = ["MeV"]
plotTitle = [" "]
plotMean = ["<"+coord+">"]
unitsY = ["mm"]
if (coord == "vertex"):
    plotXtitle = ["vertex extrapolation distance"]
    unitsX = ["m"]

#Make new canvas for plots 
w = 1800
h = 800
c = TCanvas("c", "c", w, h)
c.SetWindowSize(w + (w - c.GetWw()), h + (h - c.GetWh()))
c.Divide(2,1)
#Keep legend, histots and TFiles in scope 
legendArray=[]
histArray=[]
fileArray=[]
profileArray=[]
graphArray=[]
meanArray=[] # for the final FoM shift-nominal 

rows, cols = (len(stationName), stateN) 
profile2DArray = [[0]*cols]*rows 

i_total=0 # canvas id counter 
for i_station in range(0, len(stationName)):
    #canvas pad and legend per station 
    c.cd(i_total+1)
    #if(coord == "radial" or coord == "vertical"):
    legend =  TLegend(0.30,0.12,0.45,0.45)
    # if(coord == "vertex"):
    #     legend =  TLegend(0.10,0.45,0.45,0.9)
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
        plot.Rebin2D(2,1)
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
        profile = plot.ProfileX("profile_S"+str(i_station)+"_"+str(i_plot), 1, -1, "")
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
        tmp_holder[holder_i].GetYaxis().SetTitleOffset(1.15)
        tmp_holder[holder_i].GetYaxis().SetRangeUser(y_min, y_max) 
        tmp_holder[holder_i].GetXaxis().SetRangeUser(x_min, x_max) 
        tmp_holder[holder_i].SetTitle("S"+stationName[i_station])
        gStyle.SetTitleY(0.95)
        tmp_holder[holder_i].GetXaxis().SetTitle(plotXtitle[0]+"/ "+str(round(binN_X,3))+" "+ unitsX[0])
        tmp_holder[holder_i].GetXaxis().CenterTitle()
        tmp_holder[holder_i].SetMarkerColor( colors[i_plot] )
        tmp_holder[holder_i].SetLineColor( colors[i_plot] )
        tmp_holder[holder_i].SetMarkerStyle(20)
        tmp_holder[holder_i].SetMarkerSize(0.5)

        if (i_plot == 0):  
            if (mode == "profile"):
                profile.Draw("E")  # profile
            if (mode == "graph"):
                gr.Draw("AP")
                #gr.SetMarkerSize(2)
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
# c.SaveAs(coord+"_"+mode+"_"+"Pos_vs_mom_timeCut.C")

profile2DArray[0]=profileArray[0:9]
profile2DArray[1]=profileArray[9:]
#Calculate the sigma bands

S_array=[]

for i_station in range(0, len(stationName)):
    print(stationName[i_station])
    #print(profile2DArray[i_station])
    for i_state in range(1, stateN):
        data=profile2DArray[i_station][0]
        sim=profile2DArray[i_station][i_state]
        #loop over bins 
        bin_number_X = data.GetNbinsX()
        #print(bin_number_X)
        i_total_bins = 0
        S_mean=0
        for i_bin_x in range(0, bin_number_X):
            bin_content = data.GetBinContent(i_bin_x)
            #skip empty bins 
            if (bin_content == 0):
                continue
            #skip bins at the start (below the cut)":
            bin_center=data.GetBinCenter(i_bin_x)
            #print(bin_center)
            if (bin_center < x_min):
                continue
            #once we hit the final bin, stop:
            if (bin_center >= x_max):
                break

            #get the bin data 
            data_r = data.GetBinContent(i_bin_x)
            sim_r = sim.GetBinContent(i_bin_x)
            data_er=data.GetBinError(i_bin_x)
            sim_er=sim.GetBinError(i_bin_x)
            # print(data_r, sim_r, data_er, sim_er)
            # calculate the sigma: S = (dR)/E2, where quadrature sum of errors 
            dR = sim_r - data_r
            E2 = np.sqrt(data_er**2 + sim_er**2)
            S = dR/E2
            S_mean+=S
            i_total_bins+=1
        
        #done with a state 
        S_mean=S_mean/i_total_bins
        S_array.append(S_mean)
        #print(i_total_bins)


rows, cols = (len(stationName), stateN-1) 
S_array_2D = [[0]*cols]*rows 
S_array_2D[0]=S_array[0:8]
S_array_2D[1]=S_array[8:]    


print("For range:", x_min, "to", x_max, "MeV:")
for i_station in range(0, len(stationName)):
    print(stationName[i_station])
    for i_state in range(0, stateN-1):
        print(names[i_state+1],": ", round(S_array_2D[i_station][i_state],3))
