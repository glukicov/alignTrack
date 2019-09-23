import os
import string
import time
import decimal
import argparse, sys
import numpy as np
import matplotlib.dates as mdate
from ROOT import TH1, TH2, TFile, TCanvas, TLine, TStyle, gROOT, gStyle, TColor, TProfile, TGraphErrors
sys.path.append("/Users/gleb/")
import rootlogon as rl
rl.SetMyStyle()

GMT = 5.0/24.0 # dt conversion 

parser = argparse.ArgumentParser(description='arguments')
parser.add_argument('-f', '--fileN', help='input ROOT file')
args = parser.parse_args()

path= "MomentumSlices/vertices/station"
# path= "Extrapolation/vertices/station"
name="/h_verticalPos_vs_GPS_timeCut"

f = TFile.Open(str(args.fileN))
if f:
    print(str(args.fileN) + "is open")
else:
    print(str(args.fileN) + "not found")

stations=("12", "18")

gps_times=[[], []]
vertical=[[], []]
vertical_error = [[], []]

for i_station, station in enumerate(stations):

    print("S"+str(station))
    plot2D = f.Get(str(path+station+name))

    totalN= int(plot2D.GetEntries())
    print("Total entries", totalN)

    #find extreme x bins with data 
    min_bin_x = plot2D.FindFirstBinAbove(0, 1)
    max_bin_x = plot2D.FindLastBinAbove(0, 1)
    n_bin_x = max_bin_x - min_bin_x # number of non-zero bins
    print("Total bins with data:",n_bin_x,"min_bin_x=",min_bin_x,"max_bin_x",max_bin_x)

    #get GPS times 
    for i_bin_x in range(min_bin_x, max_bin_x):
         x_binCentre = plot2D.GetXaxis().GetBinCenter(i_bin_x)
         mpl_date = mdate.epoch2num(int(x_binCentre))
         gps_times[i_station].append(mpl_date-GMT) # epoch-float64 -> int -> MPL data as GMP Python WTF!
    
    print("Found",len(gps_times[i_station]),"GPS records")
    
    print("Start",min(mdate.num2date(gps_times[i_station])))
    print("Start",max(mdate.num2date(gps_times[i_station])))
    # print(mdate.num2date(mpl_date))

    #get <Y> times from the profile
    profile = plot2D.ProfileX("profile_S"+str(i_station), 1, -1, "")
    for i_bin_x in range(min_bin_x, max_bin_x):
         y = profile.GetBinContent(i_bin_x)
         error = profile.GetBinError(i_bin_x)
         vertical[i_station].append(y)
         vertical_error[i_station].append(error)

    print("Found",len(vertical[i_station]),"tracker records")
    print("Max",min(vertical[i_station]))
    print("Min",max(vertical[i_station]))
    



#write to file: S12, <Y_S12>, S18, <Y_S18>
f=open("OfflineData.txt", "w+")
for i_record, record in enumerate(gps_times[i_station]):
    f.write(str(gps_times[0][i_record])+", "+str(vertical[0][i_record])+", "+str(vertical[1][i_record])+"\n")
    # f.write(str(gps_times[0][i_record])+", "+str(vertical[0][i_record])+", "+str(vertical_error[0][i_record])+", "\
    #         +str(vertical[1][i_record])+", "+str(vertical_error[1][i_record])+"\n")






