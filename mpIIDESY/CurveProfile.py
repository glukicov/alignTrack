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
import numpy.polynomial.polynomial as poly
import subprocess
def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)

from ROOT import TH1F, TH2F, TF1, TCanvas, TFile, gStyle, TPaveText, TLegend, TProfile, TGraphErrors
from decimal import *
round_to = 3
getcontext().prec = round_to


def Chi2(d, s, dd, ed, es):
    chi2 = ( ( d-s+dd )**2 ) / ( ed**2 + es**2 )
    return chi2 

parser = argparse.ArgumentParser(description='arguments')
parser.add_argument('-m', '--mode', type=str)
parser.add_argument('-c', '--coord', type=str, default="radial")
parser.add_argument('--method', type=str, default=None)
parser.add_argument('-xbins', '--nX_bins', type=int, default=1)
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

method = args.method
if not (method == "sigma" or method== "chi2" or method==None):
    print("Specify 'sigma' or 'chi2' as --method=, or leave empty")
    sys.exit()

y_min = args.y_min
y_max = args.y_max
x_min = args.x_min
x_max = args.x_max
nX_bins = args.nX_bins


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
# names = ["data (aligned)", "sim default", "sim cbo"]
# states = ["run15922-15924_align", "sim_a=-0.5e-6", "sim_truth", "sim_a=0.3e-6", "sim_a=0.4e-6", "sim_a=0.5e-6", "sim_a=1e-6"]
# names = ["data (aligned)", "simulation #it{a}=-0.5#times10^{-6} mm^{-1}", "simulation #it{a}=0 mm^{-1}", "simulation #it{a}=+0.3#times10^{-6} mm^{-1}", "simulation #it{a}=+0.4#times10^{-6} mm^{-1}", "simulation #it{a}=+0.5#times10^{-6} mm^{-1}", "simulation #it{a}=+1.0#times10^{-6} mm^{-1}"]
states = ["run15922-15924_align", "sim_a=-0.5e-6", "sim_truth", "sim_a=1e-7", "sim_a=0.2e-6", "sim_a=0.3e-6", "sim_a=0.35e-6", "sim_a=0.4e-6", "sim_a=0.45e-6", "sim_a=0.5e-6", "sim_a=0.55e-6", "sim_a=0.6e-6", "sim_a=1e-6"]
names = ["data (aligned)", "sim a=-0.5e-6", "sim a=truth  ", "sim a=0.1e-6 ", "sim a=0.2e-6 ", "sim a=0.3e-6 ", "sim a=0.35e-6 ", "sim a=0.4e-6 ", "sim a=0.45e-6 ", "sim a=0.5e-6 ", "sim a=0.55e-6 ", "sim a=0.6e-6 ","sim a=1.0e-6 "]
stateN=len(states)

#Containers to store histograms in orders as the names 

colors = [1, 2, 3 ,4 ,5, 6 ,7 ,41 ,9, 49, 46, 30, 12] 
# colors = [0, 1, 2, 4 , 3, 41 ,9, 49, 46, 30, 12] 
# marker_styles= [2, 20, 4, 23, 34, 21, 33]
marker_styles=["+", "*"]
colorsStn=["purple", "orange"]
styles = [3001, 3002]
plotName = ["h_"+coord+"Pos_vs_mom_timeCut", "h_"+coord+"Pos_vs_mom"]
if (coord == "vertex"):
    plotName = ["h_radialPos_vs_vertexExtrapDist_timeCut", "h_radialPos_vs_vertexExtrapDist"]
plotYtitle = [coord+" beam position"]
plotXtitle = ["#it{p} "]
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
legendArray_=[]
histArray=[]
fileArray=[]
profileArray=[]
graphArray=[[], []]
graphArray_=[[], []]
meanArray=[] # for the final FoM shift-nominal 

rows, cols = (len(stationName), stateN) 
profile2DArray = [[0]*cols]*rows 

i_total=0 # canvas id counter 
for i_station in range(0, len(stationName)):
    #canvas pad and legend per station 
    c.cd(i_total+1)
    #if(coord == "radial" or coord == "vertical"):
    legend =  TLegend(0.25,0.12,0.40,0.45)
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

       # print("before re-binning: ", plot.GetNbinsX())
        # Rebin2D(X,Y)
        plot.Rebin2D(nX_bins,1)
        #print("after re-binning: ", plot.GetNbinsX())

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

        #print(binNumber_actual)
        gr = TGraphErrors(binNumber_actual, x, y ,ex, ey)
        gr.SetName("TGraphErrors_"+str(i_plot)+str(i_station))
        graphArray[i_station].append(gr)

        #plot depending on the mode 
        tmp_holder.append(profile)
        tmp_holder.append(gr)
        holder_i = None
        if (mode == "profile"):
            holder_i = 0 
        if (mode == "graph"):
            holder_i = 1 
        profileArray.append(tmp_holder[holder_i])
                
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

#split per station 
profile2DArray[0]=profileArray[0:int(stateN)]
profile2DArray[1]=profileArray[int(stateN):]

# #Now plot the difference of data-sim
# cD = TCanvas("cD", "cD", w, h)
# cD.SetWindowSize(w + (w - cD.GetWw()), h + (h - cD.GetWh()))
# cD.Divide(2,1)
# i_total=0 # canvas id counter 
# for i_station in range(0, len(stationName)):
#     #canvas pad and legend per station 
#     cD.cd(i_total+1)
#     #if(coord == "radial" or coord == "vertical"):
#     legend_ =  TLegend(0.13,0.12,0.40,0.45)
#     # if(coord == "vertex"):
#     #     legend =  TLegend(0.10,0.45,0.45,0.9)
#     legendArray_.append(legend_) # stroe all to keep in scope 
#     for i_state in range(1, stateN):
#         x=array( 'f', [])
#         y=array( 'f', [])
#         ex=array( 'f', [])
#         ey=array( 'f', [])
#         data=profile2DArray[i_station][0]
#         sim=profile2DArray[i_station][i_state]

#         #loop over bins 
#         bin_number_X = data.GetNbinsX()
#         #print(bin_number_X)
#         i_total_bins = 0
#         for i_bin_x in range(0, bin_number_X):
#             bin_content = data.GetBinContent(i_bin_x)
#             #skip empty bins 
#             if (bin_content == 0):
#                 continue
#             #skip bins at the start (below the cut)":
#             bin_center=data.GetBinCenter(i_bin_x)
#             #print(bin_center)
#             if (bin_center < x_min):
#                 continue
#             #once we hit the final bin, stop:
#             if (bin_center >= x_max):
#                 break

#             #get the bin data 
#             data_r = data.GetBinContent(i_bin_x)
#             sim_r = sim.GetBinContent(i_bin_x)
#             data_er=data.GetBinError(i_bin_x)*2
#             sim_er=sim.GetBinError(i_bin_x)*2
#             y.append(sim_r-data_r)
#             ey.append(np.sqrt(data_er**2+sim_er**2))
#             x.append(data.GetBinCenter(i_bin_x))
#             ex.append(data.GetBinWidth(i_bin_x)*0.5)
#             i_total_bins+=1

#         gr_ = TGraphErrors(i_total_bins, x, y ,ex, ey)
#         gr_.SetName("TGraphErrors_"+str(i_state)+str(i_station)+"_")
#         graphArray_[i_station].append(gr_)

#           #new profile/graph y-axis tools (in case tgrpah is drawn first)
#         gr_.GetYaxis().SetTitle("#Delta (data - simulation) radial beam position [mm]")            
#         gr_.GetYaxis().SetTitleSize(0.045)            
#         gr_.GetYaxis().SetLabelSize(0.040)            
#         gr_.GetXaxis().SetTitleSize(0.045)            
#         gr_.GetXaxis().SetLabelSize(0.040)            
#         gr_.GetYaxis().CenterTitle()
#         gr_.GetYaxis().SetTitleOffset(1.15)
#         gr_.GetYaxis().SetRangeUser(y_min, y_max) 
#         gr_.GetXaxis().SetRangeUser(x_min, x_max) 
#         gr_.SetTitle("S"+stationName[i_station])
#         gStyle.SetTitleY(0.95)
#         gr_.GetXaxis().SetTitle(plotXtitle[0]+"/ "+str(round(binN_X,3))+" "+ unitsX[0])
#         gr_.GetXaxis().CenterTitle()
#         gr_.SetMarkerColor( colors[i_state] )
#         gr_.SetLineColor( colors[i_state] )
#         gr_.SetMarkerStyle(marker_styles[i_state])
#         gr_.SetMarkerSize(2.0)

#         if (i_state == 1):  
#             gr_.Draw("AP")
#                 #gr.SetMarkerSize(2)
#         else:
#             gr_.Draw("P same")
            
#         #fill legend once per state      
       
#         legenValue1 = str(names[i_state])
#         legend_.AddEntry(gr_, str(legenValue1))


#     #draw legend once per canvas 
#     legend_.Draw("same")
#     i_total+=1
#     meanArray=[]
#     legend_.SetFillStyle(0)
#     legend_.SetTextSize(0.035)

# cD.Draw()
# cD.Print("dif"+coord+"_"+"Pos_vs_mom_timeCut.png")


if(method == None):
    print("No further methods to do now")
    sys.exit()

### Chi2 method ###### 
if(method == "chi2"):
    print("Computing chi2 method")

    #generate step samples 
    dd = np.arange(-5.0, 5.0, 0.05)
    stepN = len(dd)

    #Final results for all states/stations (the min point per state )
    chi2Array=[] # to store the min Chi2 
    ddArray=[] # to store the min dd 

    for i_station in range(0, len(stationName)):
        print("S"+stationName[i_station])
        for i_state in range(1, stateN):
            data=profile2DArray[i_station][0]
            sim=profile2DArray[i_station][i_state]

            #Final sum over bins
            chi2Steps=[]
            
            #loop over steps
            for i_step in range(0, stepN):

                Chi2_dd = 0 # chi2 on the current step, summing over bins 
                #print(i_step, dd[i_step])

                #loop over bins 
                bin_number_X = data.GetNbinsX()
                #print(bin_number_X)
                i_total_bins = 0
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
                    data_er=data.GetBinError(i_bin_x)*2
                    sim_er=sim.GetBinError(i_bin_x)*2
                    # print(data_r, sim_r, data_er, sim_er)
                    # calculate Chi2 as a sum from all bins Chi2(d, s, dd, ed, es):
                    Chi2_dd += Chi2(data_r, sim_r, dd[i_step], data_er, sim_er)
                    #print("Current: ", Chi2(data_r, sim_r, dd[i_step], data_er, sim_er)) 
                    #print("Total: ", Chi2_dd) 
                    i_total_bins+=1


                #done with a step
                chi2Steps.append(Chi2_dd/i_total_bins) 
            
            #done with a state 
            #print(chi2Steps)
            #print("chi2Array len", len(chi2Steps))
            #print(states[i_state])
            #print("chi2Array min", min(chi2Steps))
            chi2Array.append(min(chi2Steps))
            #print("chi2Array min ID", np.argmin(chi2Steps))
            #print("dd", dd[np.argmin(chi2Steps)])
            ddArray.append(dd[np.argmin(chi2Steps)])

            #debug sanity plot 
            #if(i_station==0 and i_state == 4):
            # f = plt.figure(figsize=(7,7))
            # plt.title("S"+stationName[i_station]+" "+states[i_state]+r" $\chi^2$: "+str(round(min(chi2Steps),3))+ r" $\delta d$: "+str(round(dd[np.argmin(chi2Steps)],3)), fontsize=12)
            # #plt.ylabel(r"$\chi^2$", fontsize=12)
            # #plt.xlabel(r"$\delta d$", fontsize=12)
            # plt.scatter(dd, chi2Steps)
            # ax0=plt.gca()
            # ax0.spines['right'].set_color('none')
            # ax0.spines['top'].set_color('none')
            # ax0.spines['left'].set_position('center')
            # ax0.spines['bottom'].set_position('center')
            # plt.show()

        #station

    # done with the Chi2 loop
    print("chi2Array", chi2Array)
    print("ddArray", ddArray)

    rows, cols = (len(stationName), stateN-1) 
    C_array_2D = [[0]*cols]*rows 
    D_array_2D = [[0]*cols]*rows 
    C_array_2D[0]=chi2Array[0:int(stateN)-1]
    C_array_2D[1]=chi2Array[int(stateN)-1:]
    D_array_2D[0]=ddArray[0:int(stateN)-1]
    D_array_2D[1]=ddArray[int(stateN)-1:]        

    print("For range:", x_min, "to", x_max, "MeV:")
    for i_station in range(0, len(stationName)):
        print(stationName[i_station])
        for i_state in range(0, stateN-1):
            print(names[i_state+1],": ", round(C_array_2D[i_station][i_state],3), round(D_array_2D[i_station][i_state],3))

    # remove low/high a
    C_array_2D[0]=C_array_2D[0][2:]
    C_array_2D[1]=C_array_2D[1][2:]
    D_array_2D[0]=D_array_2D[0][2:]
    D_array_2D[1]=D_array_2D[1][2:]

    print(C_array_2D[0])
    print(C_array_2D[1])
    print(D_array_2D[0])
    print(D_array_2D[1])


    x_ticks =  ( [0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 1.0] ,["0.1e-6 ", "0.2e-6 ", "0.3e-6 ", "0.35e-6 ", "0.4e-6 ", "0.45e-6 ", "0.5e-6 ", "0.55e-6 ", "0.6e-6 ", "1.0e-6 "])
    print(len(x_ticks[0]))
    print(len(C_array_2D[0]))

    f = plt.figure(figsize=(7,7))
    # #Plot the 0th line 
    line = [[x_ticks[0][0]-0.1,0.0], [x_ticks[0][-1]+0.1, 0.0]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
    plt.ylabel(r"$\delta r$ [mm]", fontsize=14)
    plt.xlabel(r"Curvature (a$\times10^{-6}$) [mm$^{-1}$]", fontsize=14)
    plt.xticks(fontsize=12, rotation=0) 
    plt.yticks(fontsize=12, rotation=0)
    for i_station in range(0, len(stationName)):
        a_min=(x_ticks[1][np.argmin(np.abs(D_array_2D[i_station]))])
        a_min=round(float(a_min)*1e6,2)
        plt.plot(x_ticks[0], D_array_2D[i_station], marker=marker_styles[i_station], mew=1, markersize=14, color=colorsStn[i_station], label="S"+stationName[i_station]+": $a$="+str(a_min)+r"$\times10^{-6}$", linestyle=":")

    axes=plt.gca()
    #plt.title("Range: "+str(x_min)+" to "+str(x_max)+" [MeV]; bins="+str(nX_bins), fontsize=12)
    axes.set_xlim(x_ticks[0][0]-0.1, x_ticks[0][-1]+0.1)
    axes.tick_params(axis='x',which='minor', bottom=True, top=True, direction='inout')
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
    plt.minorticks_on()
    plt.tight_layout()
    axes.legend(loc='upper center', prop={'size': 14}) # outside (R) of the plot 
    plt.savefig("D.png", dpi=250)

    f = plt.figure(figsize=(7,7))
    plt.ylabel(r"$\chi^2/\mathrm{DoF}$", fontsize=14)
    plt.xlabel(r"Curvature (a$\times10^{-6}$) [mm$^{-1}$]", fontsize=14)
    plt.xticks(fontsize=12, rotation=0) 
    plt.yticks(fontsize=12, rotation=0)
    for i_station in range(0, len(stationName)):
        a_min=(x_ticks[1][np.argmin(np.abs(C_array_2D[i_station]))])
        a_min=round(float(a_min)*1e6,2)
        plt.plot(x_ticks[0], C_array_2D[i_station], marker=marker_styles[i_station], mew=1, markersize=14, color=colorsStn[i_station], label="S"+stationName[i_station]+": $a$="+str(a_min)+r"$\times10^{-6}$", linestyle=":")

    axes=plt.gca()
    #plt.title("Range: "+str(x_min)+" to "+str(x_max)+" [MeV]; bins="+str(nX_bins), fontsize=12)
    axes.set_xlim(x_ticks[0][0]-0.1, x_ticks[0][-1]+0.1)
    axes.tick_params(axis='x',which='minor', bottom=True, top=True, direction='inout')
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
    plt.minorticks_on()
    plt.tight_layout()
    axes.legend(loc='upper center', prop={'size': 14}) # outside (R) of the plot 
    plt.savefig("Chi2.png", dpi=250)

### Sigma method ###### 
if(method == "sigma"):

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
                data_r = data.GetBinContent(i_bin_x) #if need correction: - data.GetMean(2)
                sim_r = sim.GetBinContent(i_bin_x) # if need correction :- sim.GetMean(2)
                data_er=data.GetBinError(i_bin_x)*2
                sim_er=sim.GetBinError(i_bin_x)*2
                # print(data_r, sim_r, data_er, sim_er)
                # calculate the sigma: S = (dR)/E2, where quadrature sum of errors 
                dR = sim_r - data_r
                E2 = np.sqrt(data_er**2 + sim_er**2)
                S = dR/E2
                #print(S)
                S_mean+=S
                i_total_bins+=1
            
            #done with a state 
            S_mean=S_mean/i_total_bins
            S_array.append(S_mean)
            #print(i_total_bins)

    print("Computing sigma method")
    rows, cols = (len(stationName), stateN-1) 
    S_array_2D = [[0]*cols]*rows 
    S_array_2D[0]=S_array[0:int(stateN)-1]
    S_array_2D[1]=S_array[int(stateN)-1:]    

    print("For range:", x_min, "to", x_max, "MeV:")
    for i_station in range(0, len(stationName)):
        print(stationName[i_station])
        for i_state in range(0, stateN-1):
            print(names[i_state+1],": ", round(S_array_2D[i_station][i_state],3))

    # remove low/high a
    S_array_2D[0]=S_array[2:int(stateN)-2]
    S_array_2D[1]=S_array[int(stateN)+1:-1]


    x_ticks =  ( [0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6] ,["0.1e-6 ", "0.2e-6 ", "0.3e-6 ", "0.35e-6 ", "0.4e-6 ", "0.45e-6 ", "0.5e-6 ", "0.55e-6 ", "0.6e-6 "])
    f = plt.figure(figsize=(7,7))
    #Plot the 0th line 
    line = [[x_ticks[0][0]-0.1,0.0], [x_ticks[0][-1]+0.1, 0.0]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
    plt.ylabel(r"$\sigma/\mathrm{DoF}$", fontsize=16)
    plt.xlabel(r"Curvature (a$\times10^{-6}$) [mm$^{-1}$]", fontsize=14)
    plt.xticks(fontsize=12, rotation=0) 
    plt.yticks(fontsize=12, rotation=0)
    for i_station in range(0, len(stationName)):
        x_new = np.linspace(float(min(x_ticks[0])), float(max(x_ticks[0])), num=1000) # generate x-points for evaluation 
        coefs = poly.polyfit(x_ticks[0], S_array_2D[i_station], 1) # x1 line
        ffit = poly.polyval(x_new, coefs) # plot over generated points 
        plt.plot(x_new, ffit, color=colorsStn[i_station], linestyle=":")
        x_0 = -coefs[0]/coefs[1] # x(y=0) = -c/m 
        plt.scatter(x_ticks[0], S_array_2D[i_station], marker=marker_styles[i_station], color=colorsStn[i_station], s=80, label="S"+stationName[i_station]+": $a$="+str(round(x_0,1))+r"$\times10^{-6}$")
        #print(S_array_2D[i_station])

    axes=plt.gca()
    #plt.title("Range: "+str(x_min)+" to "+str(x_max)+" [MeV]", fontsize=12)
    axes.set_xlim(x_ticks[0][0]-0.1, x_ticks[0][-1]+0.1)
    # axes.set_ylim(-3, 3)
    axes.tick_params(axis='x',which='minor', bottom=True, top=True, direction='inout')
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
    plt.minorticks_on()
    axes.legend(loc='upper left', prop={'size': 14}) # outside (R) of the plot 
    plt.tight_layout()
    plt.savefig("Sigma.png", dpi=250)