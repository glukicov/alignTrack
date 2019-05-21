# Juoyer ROOT import 
import sys
# ROOT includes 
from ROOT import TH1D, TH2D, TF1, TCanvas, TFile, gStyle, TLegend, gROOT
import numpy as np
import math
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt #for plotting  
import numpy.polynomial.polynomial as poly
import matplotlib.ticker as plticker
   
#Define constant paths and labels 
# topDir = "/Users/gleb/software/alignTrack/mpIIDESY/MomSlices/"
# topDir = ["/Users/gleb/software/alignTrack/mpIIDESY/Curve_data/", "/Users/gleb/software/alignTrack/mpIIDESY/MomSlices/"]
topDir = ["/Users/gleb/software/alignTrack/mpIIDESY/Curve_data/", "/Users/gleb/software/alignTrack/mpIIDESY/Curve_data/"]
momSlice = ("_0_800" , "_800_1000", "_1000_1500", "_1500_2300", "_2300_3100")
momName = ("0-800" , "800-1000", "1000-1500", "1500-2300", "2300-3100")
# states=("Truth", "a=+1e-6")
# statesNames=["run15922", "Truth"]
statesNames=["run15922", "run15922_Aligned"]
stateN = len(statesNames)
slicesN = len(momSlice)
fileName="gm2tracker_MomSlices_ana.root"
plotPath= "MomentumSlices/vertices/"
stations=[12, 18]
stationN=len(stations)

round_to = 3

#Final plots and canvases names (looped over i_plot)
canvasTitle = "Radial"
# plotName = ["h_radialPos_vs_time", "h_radialPos"]
plotName = ["h_radialPos_vs_time", "h_radialPos_vs_time"]
result = "<x>"
xTitle= "Radial Beam Position"

#Legen labels (looped over i_state)
legendName = statesNames
colorLine = [2, 9, 6, 5] # red, blue, green, black, purple, cyan, yellow 
colorHisto = colorLine

#Open TFiles (looped over i_state)

fileArray = [] # keep files in scope
histArray = [] # keep hists in scope
legendArray = [] #keep legends in scope 
canvasArray=[]

#Global empty containers to be filled for vertical or radial (in the main loop)
mean = -1
mean_error = -1 
sd = -1 
sd_error = -1 
meanArray=[] # for the final FoM shift-nominal 

#for plotting summary plots
meanAll = [[[0 for i_mom in range (0, slicesN) ]for i_state in range (stateN)]  for i_station in range(stationN) ] 
meanAll_error = [[[0 for i_mom in range (0, slicesN) ]for i_state in range (stateN)]  for i_station in range(stationN) ] 


###### Plotting ##########
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetLegendBorderSize(0)
# gStyle.SetLabelSize(.05, "XY")

#rad ver
canvas = TCanvas("can", "", 5830, 4100)
canvas.Divide(2,3)
i_plot = 0 
for i_mom in range(0, slicesN):
    legend =  TLegend(0.12,0.55,0.48,0.89)
    legend.SetFillStyle(0)
    i_plot+=1 
    canvas.cd(i_plot) # cd for each station and rad/ver
    i_colour=0
    for i_state in range(0, stateN): 
        for i_station in range(0, stationN):
            
            filePath = topDir[i_state] + statesNames[i_state] + "/" + fileName
            fullHistopath = plotPath + "station" +str(stations[i_station]) + "/" + str(momSlice[i_mom]) + "/" + plotName[i_state]

            # print(filePath)
            # print(fullHistopath)

            # #data
            # if(i_state==0):
            #Get the TH2F 
            rootFile = TFile.Open(filePath)
            fileArray.append(rootFile)
            histo_2D = rootFile.Get(fullHistopath)
            #Apply 30 us time cut 
            first_bin = histo_2D.GetXaxis().FindBin(30.0)
            tmpNameTH1 = "tmpNameTH1_"+str(i_mom)+str(i_station)+str(i_state) # assign a new "name pointer" to the TH1 object for each loop 
            #Get the TH1F
            hist_1D = histo_2D.ProjectionY(tmpNameTH1, first_bin, -1)
        
            # #simulation 
            # if(i_state==1):
            #     #Get the TH1F 
            #     rootFile = TFile.Open(filePath)
            #     fileArray.append(rootFile)
            #     #Get the TH1F 
            #     hist_1D = rootFile.Get(fullHistopath)


            histArray.append(hist_1D)

            #Rebin and minipulate the histo 
            if (i_plot==5):
                hist_1D.Rebin(2)
            hist_1D.GetXaxis().SetRangeUser(-50, 50) # applying a maximum range cut 
            binN=hist_1D.GetBinWidth(1)

            hist_1D.SetTitle("")
            hist_1D.GetYaxis().SetTitle("Tracks / "+str(binN)+" mm")
            hist_1D.GetXaxis().SetTitle(xTitle +" [mm]")       
            hist_1D.GetYaxis().SetTitleOffset(1.4)
            hist_1D.GetXaxis().SetTitleOffset(1.4)
            hist_1D.GetYaxis().CenterTitle()
            hist_1D.GetXaxis().CenterTitle()
            hist_1D.GetXaxis().SetTitleSize(0.04)
            hist_1D.GetYaxis().SetTitleSize(0.04)

            #Draw on canvas 
            hist_1D.SetLineColor(colorHisto[i_colour])
            hist_1D.SetMarkerSize(4)
    #                 hist_1D.SetLineSize(4)
            if (i_colour == 0):
                hist_1D.Draw("E1")

            else:
                hist_1D.Draw("E1 same")

            #Get stats from hist 
            mean = hist_1D.GetMean()
            # print(statesNames[i_state])
            # print(stations[i_station])
            # print(momName[i_mom])
            # print(mean)
            
            mean_error = hist_1D.GetMeanError()
            sd = hist_1D.GetRMS()
            sd_error = hist_1D.GetRMSError()
            legenObject = hist_1D
            meanArray.append(mean)
            meanAll[i_station][i_state][i_mom]=(mean)
            meanAll_error[i_station][i_state][i_mom]=(mean_error)

            #take care of the legend
            legenValue1 = str("S"+str(stations[i_station])+" "+legendName[i_state])+": "+str(result)+": "+str(round(mean,round_to))+" #pm "+str(round(mean_error, round_to)) 
            legenValue2 = "#sigma: "+str(round(sd,round_to))+" #pm "+str(round(sd_error,round_to))
            legend.AddEntry(legenObject,"#splitline{"+str(legenValue1)+"}{           "+str(legenValue2)+"}","L") # make appropriate spacing 
            legend.SetTextSize(.038)
            legend.Draw("same")

            i_colour=i_colour+1
        #Do some final massaging per pad 
        legend.SetHeader(str(canvasTitle)+": P=" + momName[i_mom] + " MeV", "C"); # option "C" allows to center the header
        legendArray.append(legend)
        meanArray=[]


#one canvas per rad/ver station 
canvas.Draw()
canvas.Print("Extrap.png")

#shift curves to 0
# for i_station in range(0, stationN):
#     for i_state in range(1, 2):
#         for i_mom in range(0, slicesN):
#             meanAll[i_station][i_state][i_mom]=meanAll[i_station][i_state][i_mom]+5.5
# print(meanAll)


xLabel="<P> in range [MeV]"
yLabel =r"$\Delta$ Mean Radial [mm]"
x_ticks = (400, 900, 1250, 1900, 2700)
aExtra = (0, 0, 0, 0, -200)
a_labels = momName
colors2D = [["blue", "purple"], ["red", "orange"]]
labels = ["S12", "S18"]
# labels2D = [["S12 run 15922", "S12 truth"],["S18 run 15922", "S18 truth"]]
labels2D = [["S12 run 15922", "S12 run 15922 aligned"],["S18 run 15922", "S18 run 15922 aligned"]]
i_plot=0
fig = plt.figure(figsize=(19,14) )
plt.ylim(2, 11)
for i_state in range(0, 2):
        axes = plt.gca()
        axes.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.ylabel(yLabel, fontsize=24)
        plt.xlabel(xLabel, fontsize=24)
        plt.xticks(fontsize=24, rotation=0) 
        plt.yticks(fontsize=24, rotation=0)
        plt.minorticks_on()
        axes.tick_params(axis='x', which='minor',bottom=False)
        axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
        #Plot data
        for i_station in range(0, stationN):
            plt.plot(x_ticks,  meanAll[i_station][i_state], color=colors2D[i_station][i_state], marker=".", linewidth=0)  
            plt.errorbar(x_ticks ,meanAll[i_station][i_state],  yerr=meanAll_error[i_station][i_state], color=colors2D[i_station][i_state], label=labels2D[i_station][i_state], elinewidth=2, linewidth=0)  
            #Fit a line for dR only 
            # print("Plotting curve", i_station, i_state, colors2D[i_station][i_state])
            x_new = np.linspace(float(min(x_ticks)), float(max(x_ticks)), num=1000) # generate x-points for evaluation 
            coefs = poly.polyfit(x_ticks, meanAll[i_station][i_state], 2) # x2 curve
            ffit = poly.polyval(x_new, coefs) # plot over generated points 
            plt.plot(x_new, ffit, color=colors2D[i_station][i_state], linewidth=2)
           

#anotate points (once per station)
for i_point, txt in enumerate(a_labels):
    axes.annotate(txt, (x_ticks[i_point]+aExtra[i_point], meanAll[i_station][i_state][i_point]), fontsize=20)
    
axes.legend(loc='lower center', fontsize=28)
    

plt.tight_layout()
plt.savefig("Summary_Extrap.png", dpi=500)