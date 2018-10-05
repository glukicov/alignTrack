#!/usr/bin/env python

#Plotter

from ROOT import *

Tf = TFile('debug.root', 'RECREATE')
c = TCanvas("c", "HitsMom", 200, 10, 700, 500)
c.Divide(1,1)
HitsMom= TH2F("HitsMom", ";number of hits; Momentum [MeV]",  32, 0.0, 32.0, 35, 0.0, 3500.0) 

#Read file and fill histos 
with open("debug.txt") as f:
    for line in f:  #Line is a string
        number = line.split()
        if (number[0] != "0") :
       		HitsMom.Fill(float(number[0]), float(number[1])) 

c.cd(1)
HitsMom.Draw("COLZ")

c.Print("HitsMom.png")
c.Print("HitsMom.root")


