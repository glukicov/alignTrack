# Juoyer ROOT import 
import sys
sys.path.append("/usr/local/Cellar/root/6.14.04_2/lib/root")
# ROOT includes 
from ROOT import TH1F, TH2F, TF1, TCanvas, TFile, gStyle, TLegend, gROOT

#Define constant paths and labels 
station12Path = "Extrapolation/vertices/station12/"
station18Path = "Extrapolation/vertices/station18/"
scr_nominal = "/Users/gleb/software/alignTrack/mpIIDESY/15922/trackRecoPlots_15922_Iter1.root"
scr_shift = "/Users/gleb/software/alignTrack/mpIIDESY/15922_Iter2/trackRecoPlots_15922_Iter2.root"
stateN = 2 # nominal (0) and shift (1)

#Final plots and canvases names (looped over i_plot)
canvasName = ["s12_vertical", "s18_vertical", "s12_radial", "s18_radial"]
canvasTitle = ["S12 Vertical", "S18 Vertical", "S12 Radial", "S18 Radial"]
plotPath = [station12Path, station18Path, station12Path, station18Path]
plotName=["h_verticalPos_vs_time", "h_verticalPos_vs_time", "h_radialPos_vs_time", "h_radialPos_vs_time"] # TFile TH2 names  

#Legen labels (looped over i_state)
legendName = ["Iteration1", "Iteration2"]
colorLine = [1, 8] # black, green
colorHisto = [2, 5] # red, yellow 

#Open TFiles (looped over i_state)
scrFile_nominal = TFile.Open(scr_nominal)
scrFile_shift = TFile.Open(scr_shift)
fileName = [scrFile_nominal, scrFile_shift] # nominal (0) and shift (1)

#Global empty containers to be filled for vertical or radial (in the main loop)
result = "-1"
mean = -1
mean_error = -1 
sd = -1 
sd_error = -1 
histo1D = []

###### Plotting ##########

for i_plot in range(0, len(canvasName)):

    #Make new canvas for plots (4 in total)
    canvas = TCanvas(str(canvasName[i_plot]), " ", 800, 600)
    canvas.Divide(1)
    canvas.cd(1)
    
    #Postion legend once per main loop accordingly 
    if( canvasName[i_plot] == "s12_vertical" or  canvasName[i_plot] == "s18_vertical"):
        legend =  TLegend(0.87, 0.87, 0.59, 0.55)
    if( canvasName[i_plot] == "s12_radial" or  canvasName[i_plot] == "s18_radial"):
        legend =  TLegend(0.105, 0.89, 0.53, 0.55)   
           
    #Loop over the 2 iterations (nominal and shift)
    for i_state in range(0, stateN):
    
        #Get the TH2F 
        histo_2D = fileName[i_state].Get(str(plotPath[i_plot])+plotName[i_plot])
        
        #Apply 30 us time cut 
        first_bin = histo_2D.GetXaxis().FindBin(30.0)
        
        # Project the beam postion into 1D histo
        hist_1D = histo_2D.ProjectionY("", first_bin, -1)
        hist_1D.SetTitle("")
        hist_1D.GetXaxis().SetRangeUser(-70, 70) # applying a maximum range cut 
        norm = hist_1D.GetEntries();
        hist_1D.Scale(1/norm); # normalise the histo 
        hist_1D.GetYaxis().SetTitle("Enteries (Normalised 1/N)")
        hist_1D.GetYaxis().SetTitleOffset(1.4);
            
        #Do a fit for vertical only 
        if( canvasName[i_plot] == "s12_vertical" or  canvasName[i_plot] == "s18_vertical"):
            
            #Define a Gaussian fit function between -30 to 30 mm 
            gF = TF1("gF", "gaus", -30.0, 30.0)
            hist_1D.Fit(gF, "QR", ) # quite fit over the specified range above
            gF.SetLineColor(colorLine[i_state])
            gF.Draw("same")

            mean = round(gF.GetParameter(1), 2)
            mean_error = round(gF.GetParError(1), 2)
            sd = round(gF.GetParameter(2), 2)
            sd_error = round(gF.GetParError(2), 2)
            result = "<y>"
            legenObject = gF

        #Extract mean from a radial histogram only  
        if( canvasName[i_plot] == "s12_radial" or  canvasName[i_plot] == "s18_radial"):
            
            #hist_1D.SetFillColor(colorIndex[i_state])
            hist_1D.Draw("E1")
            hist_1D.Draw("same")
            mean = round(hist_1D.GetMean(), 2)
            mean_error = round(hist_1D.GetMeanError(), 2)
            sd = round(hist_1D.GetRMS(), 2)
            sd_error = round(hist_1D.GetRMSError(), 2)
            result = "<x>"
            legenObject = hist_1D
            
            
        #take care of the legend
        legend.SetHeader(canvasTitle[i_plot], "C"); # option "C" allows to center the header
        legenValue1 = str(legendName[i_state])+": "+str(result)+": "+str(mean)+" #pm "+str(mean_error) 
        legenValue2 = "#sigma: "+str(sd)+" #pm "+str(sd_error)
        legend.AddEntry(legenObject,"#splitline{"+str(legenValue1)+"}{                 "+str(legenValue2)+"}","L") # make appropriate spacing 
        legend.SetTextSize(.028)
        legend.Draw("same")

    #Do some final massagin and save to a file 
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.023)
    gROOT.ForceStyle()
    canvas.Draw()
    canvas.Print(str(canvasName[i_plot])+".png")