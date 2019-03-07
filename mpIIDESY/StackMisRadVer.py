#!/usr/bin/python

####################################################################
# Plot Radial and Vertical Misalignment given the FHICL file 
#####################################################################
import argparse, sys, os

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-path')
parser.add_argument('-label')
parser.add_argument('-mode')
parser.add_argument('-scale')
parser.add_argument("-cases", nargs='+')
args = parser.parse_args()

import matplotlib.pyplot as plt #for plotting 
import matplotlib.ticker as ticker
import numpy as np  # smart arrays 
import itertools # smart lines
import subprocess
import pandas as pd
import re
from ROOT import TCanvas, TH1, TLegend, TFile, gStyle

# CONSTANTS 
moduleN = 8 #number of movable detectors
stationN = 3 # S0, S12, S18 
globalN = 2 # X, Y
parN = stationN * moduleN

def getOffsets(f, name):
	offsets = [] #tmp storage buffer 
	for line in f:
		if re.match(name, line):
			copy = True
			offsets=line 

		else:
			copy = False

	return offsets   
		

cases=args.cases
path = str(args.path)
label=str(args.label)
mode=str(args.mode)
scale=str(args.scale)

colours = ["green", "blue", "black", "orange", "purple", "red", "brown", "grey"]
markers = [".", "v", "P", "x", "D", "^", "2", "p"]
moduleArray=range(1, moduleN+1)

totalR=[]
totalH=[]
if (mode == "mis"):
	for i_total, i_case in enumerate(cases):
		fullPath=path+"/"+str(i_case)
		file=fullPath+"/RunTrackingDAQ_align.fcl"
		#print(fullPath) 

			###Store the offsets from file
		f=open(file, "r")
		radial = (getOffsets(f, "services.Geometry.strawtracker.strawModuleRShift12:"))
		radial = radial.replace("services.Geometry.strawtracker.strawModuleRShift12: [", " ") 
		radial = radial.replace("]", "") 
		radialOff = np.array([float(r) for r in radial.split(',')])
		radialOff=radialOff[0:8] 
		radialOff=radialOff*1e3 # mm -> um 
		totalR.append(radialOff)
		#print("Radial:", radialOff)

		f=open(file, "r")
		vertical = (getOffsets(f, "services.Geometry.strawtracker.strawModuleHShift12:"))
		vertical = vertical .replace("services.Geometry.strawtracker.strawModuleHShift12: [", "") 
		vertical = vertical .replace("]", "") 
		verticalOff = np.array([float(v) for v in vertical .split(',')])
		verticalOff=verticalOff[0:8] 
		verticalOff=verticalOff*1e3 # mm -> um 
		totalH.append(verticalOff)
		#print("Vertical:",  verticalOff)



	yMin = -200
	yMax = 200
	plt.subplot(211) # X 
	plt.rcParams.update({'font.size': 10})
	axes = plt.gca()
	axes.set_xlim(0.5, 8.5)
	axes.set_ylim(yMin, yMax)
	plt.title("Misalignment: Radial "+str(label), fontsize=10)
	plt.ylabel("Misalignment [um]", fontsize=10)
	line = [[0.5,0.0], [8.5, 0.0]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	for i_total, i_case in enumerate(cases):
		plt.plot(moduleArray, totalR[i_total], marker=markers[i_total], color=colours[i_total], label=str(i_case))
	for i_module in range(0, 8):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')

	plt.xlabel("Module", fontsize=10)
	axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	# plt.legend(loc='center')

	plt.subplot(212) # Y 
	axes = plt.gca()
	axes.set_xlim(0.5, 8.5)
	axes.set_ylim(yMin, yMax)
	plt.title("Misalignment: Vertical "+str(label), fontsize=10)
	plt.ylabel("Misalignment [um]", fontsize=10)
	line = [[0.5,0.0], [8.5, 0.0]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	for i_total, i_case in enumerate(cases):
		plt.plot(moduleArray, totalH[i_total], marker=markers[i_total], color=colours[i_total], label=str(i_case))
	for i_module in range(0, 8):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')


	plt.xlabel("Module", fontsize=10)
	plt.tight_layout()
	plt.subplots_adjust(right=0.8)
	axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.savefig("Mis_"+str(label)+".png", dpi=600)

if (mode == "plot"):
	#hack in the "0th" case at the front 
	cases.insert(0, 0) 

	station12Path = "Extrapolation/vertices/station12/"
	hitTimeHistoPath_s12 = "HitSummary/Station_12/h_hitTime"

	canvasName = "s12_vertical"
	canvasTitle = "S12"
	plotPath = station12Path
	plotName="h_verticalPos_vs_time"
	legendName = cases 
	hitsTotal_s12_nom = - 1 # global container 

	#create main canvas 
	canvas = TCanvas(canvasName, canvasTitle, 16000, 10000)
	canvas.Divide(2,2)
	legend =  TLegend(0.87, 0.87, 0.59, 0.55) 
	legend_tracks =  TLegend(0.87, 0.87, 0.59, 0.55) 
	legend_hits =  TLegend(0.87, 0.87, 0.59, 0.55) 
	
	scrFiles = [] # keep all opened files in scope 
	clonedHistos=[]

	#now stack the dist.
	for i_total, i_case in enumerate(cases):
		
		#open the case root file 
		fullPath=path+"/"+str(i_case)
		fileName=fullPath+"/trackRecoPlots.root"
		scrFile = TFile.Open(fileName)
		scrFiles.append(scrFile) # !!! need to have ROOT files open/stored in memory otherwise TH1 objects go "out-of-scope"
			   
		#Get the TH2F 
		histo_2D = scrFile.Get(plotPath+plotName)
		 
		#Apply 30 us time cut and get TH1
		first_bin = histo_2D.GetXaxis().FindBin(30.0)
		tmpNameTH1 = "tmpNameTH1_"+str(i_case)+str(i_total) # assign a new "name pointer" to the TH1 object for each loop 
		hist_1D = histo_2D.ProjectionY(tmpNameTH1, first_bin, -1)
		hist_1D.Rebin(2)
		hist_1D.SetTitle("")
		hist_1D.GetXaxis().SetRangeUser(-50, 50) # applying a maximum range cut 
		hist_1D.GetYaxis().SetTitleOffset(1.4);
		hist_1D.GetXaxis().SetTitleOffset(1.4);
		hist_1D.SetLineColor(i_total+1) # 0=white... 
		
		#Normalisation 
		#only for nominal case
		if (i_case == 0):
			hits_histo_nom = scrFile.Get(hitTimeHistoPath_s12)
			first_bin_nom = hits_histo_nom.GetXaxis().FindBin(0.03) # get integral for 30us ->
			last_bin_nom = hits_histo_nom.GetXaxis().FindBin(1.0) # get integral for 30us ->  
			hitsTotal_s12_nom=hits_histo_nom.Integral(first_bin_nom, last_bin_nom) 

		#all cases (inc. nominal)
		hits_histo = scrFile.Get(hitTimeHistoPath_s12)
		first_bin = hits_histo.GetXaxis().FindBin(0.03) # get integral for 30us ->
		last_bin = hits_histo.GetXaxis().FindBin(1.0) # get integral for 30us ->  
		hitsTotal_s12=hits_histo.Integral(first_bin, last_bin) 

		tracks = hist_1D.GetEntries()
		hits = hitsTotal_s12/hitsTotal_s12_nom
		# norm =  tracks/hits	  
		norm =  hits	  

		#unscaled  
		#if (scale == "none"):
		canvas.cd(1)
		legend.SetHeader(canvasTitle+" "+str(label)+" [unnorm.] cases: ", "C")
		hist_1D.GetYaxis().SetTitle("Enteries")
		legenObject = hist_1D	
		legenValue1 = str(legendName[i_total]) + " <y>= " + str(round(hist_1D.GetMean(), 3))
		legend.AddEntry(legenObject,str(legenValue1),"L") 
		legend.SetTextSize(.028)
		if (i_case == 0):
			hist_1D.Draw("")
			legend.Draw("")
			
		else:
			hist_1D.Draw("same")
			legend.Draw("same")

		#scaled by tracks 
		# #if (scale == "tracks"):
		canvas.cd(2)
		legend_tracks.SetHeader(canvasTitle+" "+str(label)+" [norm. tracks] cases: ", "C")
		cloneNameTH1 = "hist_1D_tracks_"+str(i_case)+str(i_total)
		hist_1D_tracks = hist_1D.Clone(cloneNameTH1)
		clonedHistos.append(hist_1D_tracks)
		hist_1D_tracks.Scale(1/tracks); # normalise the histo 
		hist_1D_tracks.GetYaxis().SetTitle("Enteries (Normalised: 1/Tracks)")
		legenObject = hist_1D_tracks	
		legenValue1 = str(legendName[i_total]) + " <y>= " + str(round(hist_1D_tracks.GetMean(), 3))
		legend_tracks.AddEntry(legenObject,str(legenValue1),"L") 
		legend_tracks.SetTextSize(.028)
		if (i_case == 0):
			hist_1D_tracks.Draw("")
			legend_tracks.Draw("")
			
		else:
			hist_1D_tracks.Draw("same")
			legend_tracks.Draw("same")

		#scaled by (case_hits/nominal_hits)/Tracks 
		#if (scale == "hits"):
		canvas.cd(3)
		legend_hits.SetHeader(canvasTitle+" "+str(label)+" [norm. hits] cases: ", "C")
		cloneNameTH1 = "hist_1D_hits_"+str(i_case)+str(i_total)
		hist_1D_hits = hist_1D.Clone(cloneNameTH1)
		clonedHistos.append(hist_1D_hits)
		hist_1D_hits.Scale(1/norm); # normalise the histo 
		hist_1D_hits.GetYaxis().SetTitle("Enteries (Normalised: 1/(case_hits/nominal_hits))")
		legenObject = hist_1D_hits	
		legenValue1 = str(legendName[i_total]) + " <y>= " + str(round(hist_1D_hits.GetMean(), 3))
		legend_hits.AddEntry(legenObject,str(legenValue1),"L") 
		legend_hits.SetTextSize(.028)
		if (i_case == 0):
			hist_1D_hits.Draw("")
			legend_hits.Draw("")
			
		else:
			hist_1D_hits.Draw("same")
			legend_hits.Draw("same")

			
	#Do some final massagin and save to a file
	gStyle.SetOptStat(0)
	gStyle.SetOptFit(0)
	gStyle.SetLegendBorderSize(0)
	gStyle.SetLegendTextSize(0.023)
	canvas.Draw()
	canvas.Print(str(canvasName)+str(label)+".png")
