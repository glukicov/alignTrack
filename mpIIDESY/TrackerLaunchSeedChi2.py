#!/usr/local/bin/python2.7

####################################################################
# FoM Plots for comparison of actual misalignment vs PEDE results 
#
# 
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################



from ROOT import *



Chi2Recon=[]
Chi2Est=[]

f = open("chi2.txt")
for line in f:  #Line is a string
    number_str = line.split()
    Chi2Recon.append(float(number_str[0]))
    Chi2Est.append(float(number_str[1]))
    print line
        

# with open("chi2.txt", "r") as ins:
#     for line in ins:
#         Chi2Recon.append(float(line[0]))
#         Chi2Est.append(float(line[1]))
#         print line[0]
     


##################PLOTING##############################

f = TFile('Chi2Tracker.root','RECREATE')

h_Chi2Recon  = TH1F("h_Chi2Recon", "Recon Chi 2", 119, 151, 154)
h_Chi2Est  = TH1F("h_Chi2Est", "Est Chi 2", 39, 150.1, 150.2)

for i in range(0, len(Chi2Recon)):
	h_Chi2Recon.Fill(Chi2Recon[i])
	h_Chi2Est.Fill(Chi2Est[i])



f.Write()
f.Close()

print "Done!"





		
