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
topDir = "/Users/gleb/software/alignTrack/mpIIDESY/AlignCurvature/"
momSlice = ( "_0_737", "_737_862", "_862_937", "_937_1037", "_1037_1112", "_1112_1212", "_1212_1287", "_1287_1362", "_1362_1462", "_1462_1587", "_1587_1687", "_1687_1862", "_1862_2087", "_2087_2362", "_2362_3100")
statesNames=["run15922_aligned_noVolumes", "sim_truth", "sim_a=0.5e-6", "sim_a=0.4e-6", "sim_a=0.3e-6", "sim_a=0.2e-6", "sim_a=1e-7", "sim_a=-0.5e-6"]
stateN = len(statesNames)
slicesN = len(momSlice)
fileName="gm2tracker_MomSlices_ana.root"
plotPath= "MomentumSlices/vertices/"
stations=[12, 18]
stationN=len(stations)


#Final plots and canvases names (looped over i_plot)
plotName = ["h_radialPos_vs_time", "h_radialPos"] 




#for plotting summary plots
meanAll = [[[0 for i_mom in range (0, slicesN) ]for i_state in range (stateN)]  for i_station in range(stationN) ] 
meanAll_error = [[[0 for i_mom in range (0, slicesN) ]for i_state in range (stateN)]  for i_station in range(stationN) ] 

data_shift = [[0 for i_state in range (stateN)]  for i_station in range(stationN) ]
###### Read-in data ##########
for i_mom in range(0, slicesN):
    for i_state in range(0, stateN): 
        for i_station in range(0, stationN):
            
            filePath = topDir + statesNames[i_state] + "/" + fileName
            
            #need time-cut for data 
            if (i_state == 0):
                fullHistopath = plotPath + "station" +str(stations[i_station]) + "/" + str(momSlice[i_mom]) + "/" + plotName[0]
                #Get the TH2F 
                rootFile = TFile.Open(filePath)
                histo_2D = rootFile.Get(fullHistopath)
                #Apply 30 us time cut 
                first_bin = histo_2D.GetXaxis().FindBin(30.0)
                tmpNameTH1 = "tmpNameTH1_"+str(i_mom)+str(i_station)+str(i_state) # assign a new "name pointer" to the TH1 object for each loop 
                #Get the TH1F
                hist_1D = histo_2D.ProjectionY(tmpNameTH1, first_bin, -1)

                #also get the shift for data just once from the total dist. 
                if (i_mom == 0):
                    fullHistopath = plotPath + "station" +str(stations[i_station]) + "/" + plotName[0]
                    rootFile = TFile.Open(filePath)
                    histo_2D = rootFile.Get(fullHistopath)
                    #Apply 30 us time cut 
                    first_bin = histo_2D.GetXaxis().FindBin(30.0)
                    tmpNameTH1 = "tmpNameTH1_"+str(i_mom)+str(i_station)+str(i_state) # assign a new "name pointer" to the TH1 object for each loop 
                    #Get the TH1F
                    hist_1D = histo_2D.ProjectionY(tmpNameTH1, first_bin, -1)
                    data_shift[i_station][i_state] = hist_1D.GetMean()
        
            #simulation 
            else:
                fullHistopath = plotPath + "station" +str(stations[i_station]) + "/" + str(momSlice[i_mom]) + "/" + plotName[1]
                #Get the TH1F 
                rootFile = TFile.Open(filePath)
                hist_1D = rootFile.Get(fullHistopath)

                #also get the shift for simulation just once from the total dist. 
                if (i_mom == 0):
                    fullHistopath = plotPath + "station" +str(stations[i_station]) + "/" + plotName[1]
                    #Get the TH1F 
                    rootFile = TFile.Open(filePath)
                    hist_1D = rootFile.Get(fullHistopath)
                    data_shift[i_station][i_state] = hist_1D.GetMean()


            #set range for collimator 
            hist_1D.GetXaxis().SetRangeUser(-50, 50) # applying a maximum range cut 
            binN=hist_1D.GetBinWidth(1)

            #Get stats from hist 
            mean = hist_1D.GetMean()
            mean_error = hist_1D.GetMeanError()
            meanAll[i_station][i_state][i_mom]=mean
            meanAll_error[i_station][i_state][i_mom]=mean_error


### Make summary plot ###
xLabel="<Track momentum in range> [MeV]"
yLabel ="<Radial Position> (arbitrary) [mm]"
x_ticks = (549, 799, 899, 987, 1074, 1162, 1249, 1324, 1412, 1524, 1637, 1774, 1974, 2224, 2731)
labels = ["S12", "S18"] 
colors = ["black", "red", "green", "blue", "purple", "orange", "brown", "grey"] # per state
lineStyle= ["-", ":", ":", ":", ":", ":", ":", "-."] # per state 
global_fontsize=8

fig = plt.figure(1)
for i_station in range(0, stationN):
    plt.subplot(int( "2"+"1"+str(int(i_station+1)) )) 
    axes = plt.gca()
    plt.tight_layout()
    for i_state in range(0, stateN):
        axes.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.ylabel(yLabel, fontsize=global_fontsize)
        plt.xlabel(xLabel, fontsize=global_fontsize)
        plt.xticks(fontsize=global_fontsize, rotation=0) 
        plt.yticks(fontsize=global_fontsize, rotation=0)
        plt.minorticks_on()
        axes.tick_params(axis='x', which='minor',bottom=False)
        axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')

        # shift data and simulation:
        meanAll[i_station][i_state] = np.array(meanAll[i_station][i_state]) - data_shift[i_station][i_state]

        #Plot data
        plt.plot(x_ticks,  meanAll[i_station][i_state], color=colors[i_state], marker=".", linewidth=0)  
        plt.errorbar(x_ticks ,meanAll[i_station][i_state],  yerr=meanAll_error[i_station][i_state], color=colors[i_state], label=labels[i_station]+" "+statesNames[i_state], elinewidth=2, linewidth=0)  
        #Fit a line for dR only 
        # print("Plotting curve", i_station, i_state, colors2D[i_station][i_state])
        x_new = np.linspace(float(min(x_ticks)), float(max(x_ticks)), num=1000) # generate x-points for evaluation 
        coefs = poly.polyfit(x_ticks, meanAll[i_station][i_state], 2) # x2 curve
        ffit = poly.polyval(x_new, coefs) # plot over generated points 
        plt.plot(x_new, ffit, color=colors[i_state], linewidth=2, linestyle=lineStyle[i_state])
    
    #make legend per subplot            
    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
plt.savefig("Summary_Extrap.png", dpi=500)

#data only 
colors = ["red", "blue"] # per state
fig = plt.figure(2)
for i_station in range(0, stationN):
    axes = plt.gca()
    plt.tight_layout()
    axes.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylabel(yLabel, fontsize=global_fontsize)
    plt.xlabel(xLabel, fontsize=global_fontsize)
    plt.xticks(fontsize=global_fontsize, rotation=0) 
    plt.yticks(fontsize=global_fontsize, rotation=0)
    plt.minorticks_on()
    axes.tick_params(axis='x', which='minor',bottom=False)
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')

    # shift data back:
    meanAll[i_station][0] = np.array(meanAll[i_station][0]) + data_shift[i_station][0]

    #Plot data
    plt.plot(x_ticks,  meanAll[i_station][0], color=colors[i_station], marker=".", linewidth=0)  
    plt.errorbar(x_ticks ,meanAll[i_station][0],  yerr=meanAll_error[i_station][0], color=colors[i_station], label=labels[i_station]+" "+statesNames[0], elinewidth=2, linewidth=0)  
    #Fit a line for dR only 
    # print("Plotting curve", i_station, i_state, colors2D[i_station][i_state])
    x_new = np.linspace(float(min(x_ticks)), float(max(x_ticks)), num=1000) # generate x-points for evaluation 
    coefs = poly.polyfit(x_ticks, meanAll[i_station][0], 2) # x2 curve
    ffit = poly.polyval(x_new, coefs) # plot over generated points 
    #plt.plot(x_new, ffit, color=colors[i_station], linewidth=2, linestyle=lineStyle[0])
    
    #make legend per subplot            
    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
plt.savefig("Summary_Extrap_data.png", dpi=500)