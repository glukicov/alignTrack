import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
from matplotlib.ticker import MaxNLocator
import argparse 
import collections, itertools
import sys, os
import subprocess
import datetime
import pickle
import scipy.interpolate
import scipy.ndimage
# from ROOT import TH3D, TFile, TCanvas

#Pass some commands 
arg_parser = argparse.ArgumentParser(description='Input data files')
arg_parser.add_argument('--dbFile', type=str, required=True, dest='dbFile')
arg_parser.add_argument('--gpsFile', type=str, required=True, dest='gpsFile')
arg_parser.add_argument('--bins', type=int, required=True, dest='bins')
args = arg_parser.parse_args()
dbFile = args.dbFile
gpsFile = args.gpsFile
bins = args.bins

### Constants 
db_type = [('times',np.float64),('qhv',np.float64)]
gps_type = [('times',np.float64),('s12',np.float64), ('s18',np.float64)]

# Data types 
db_collection_shape = ['time', 'qhv']
gps_collection_shape= ['time', 's12', 's18']
db_collection = collections.namedtuple('db_collection', db_collection_shape)
db_collection_clean = collections.namedtuple('db_collection_clean', db_collection_shape)
db_collection_reshaped = collections.namedtuple('db_collection_clean', db_collection_shape)
gps_collection = collections.namedtuple('gps_collection', gps_collection_shape)
gps_collection_reshaped = collections.namedtuple('gps_collection', gps_collection_shape)

stations=("S12", "S18")

def main():

    print("Starting plots on:", datetime.datetime.now())

    #Get data 
    print("Getting data from files...")
    db_data = getDBData(dbFile)
    gps_data = getGPSData(gpsFile)

    #Clean DB data (same start/end points)
    print("Cleaning DB data...")
    db_data_clean = cleanDB(db_data, gps_data)

    #Rebin data 
    print("Re-binning data into",bins,"bins")
    db_data_reshaped = dbReshape(db_data_clean, nbins=bins)
    gps_data_reshaped = gpsReshape(gps_data, nbins=bins)

    #Plot data 
    print("Plotting data..")
    plotData(db_data_reshaped, gps_data_reshaped)

    print("Finishing plots on:", datetime.datetime.now())


def getDBData(file_name):
    db_data = np.genfromtxt(file_name, dtype=db_type, delimiter=',')
    times = db_data['times']
    qhv = db_data['qhv'] 
    data = db_collection(times, qhv)
    return data

def getGPSData(file_name):
    gps_data = np.genfromtxt(file_name, dtype=gps_type, delimiter=',')
    times = gps_data['times']
    s12 = gps_data['s12'] 
    s18 = gps_data['s18'] 
    data = gps_collection(times, s12, s18)
    return data

def cleanDB(db_collection, gps_collection):
    time_db = db_collection.time
    qhv = db_collection.qhv
    time_gps = gps_collection.time

    new_time_db = [] 
    new_qhv=[]

    # remove extreme db times 
    for i_entry in range(len(time_db)):
        if(time_db[i_entry] >= time_gps[0] and time_db[i_entry] <= time_gps[-1]):
            new_time_db.append(time_db[i_entry])
            new_qhv.append(qhv[i_entry])

     
    print("DB length", len(new_time_db))
    print("GPS length", len(time_gps))
    print( "DB min, max", np.min(mdate.num2date(new_time_db)),  np.max(mdate.num2date(new_time_db)) )
    print("GPS min, max", np.min(mdate.num2date(time_gps)), np.max(mdate.num2date(time_gps)) )
    clean_data = db_collection_clean(new_time_db, new_qhv)
    return clean_data


def dbReshape(db_collection_clean, nbins=20): 
    data_len = len(db_collection_clean[0])
    rescale = int(data_len/nbins) # elements per bin 
    total_kept = rescale * nbins # can fit in bins equally
    total_removed = data_len - total_kept # drop elements that can fit 
    rescaled_data_array =[] 
    for array in db_collection_clean:  #lop over each array in the collection 
        array = (array[:-total_removed or None]) # remove elements 
        rescaled_array = [] 
        for i_bin in range(nbins):  #loop over bins, as a range of rescaled elements 
            lower = rescale*i_bin 
            upper = rescale*(i_bin+1)
            rescaled_array.append(np.mean(array[lower:upper])) # take the mean of the slice (bin)
        rescaled_data_array.append(rescaled_array)
    dataR = db_collection_reshaped(*rescaled_data_array)  # push arrays into the collection 
    return dataR

def gpsReshape(gps_collection, nbins=20):
    data_len = len(gps_collection[0])
    rescale = int(data_len/nbins) # elements per bin 
    total_kept = rescale * nbins # can fit in bins equally
    total_removed = data_len - total_kept # drop elements that can fit 
    rescaled_data_array =[] 
    for array in gps_collection:  #lop over each array in the collection 
        array = (array[:-total_removed or None]) # remove elements 
        rescaled_array = [] 
        for i_bin in range(nbins):  #loop over bins, as a range of rescaled elements 
            lower = rescale*i_bin 
            upper = rescale*(i_bin+1)
            rescaled_array.append(np.mean(array[lower:upper])) # take the mean of the slice (bin)
        rescaled_data_array.append(rescaled_array)
    dataR = gps_collection_reshaped(*rescaled_data_array)  # push arrays into the collection 
    return dataR
    

def plotData(db_collection, gps_collection):
    time_db = db_collection.time
    qhv = db_collection.qhv
    time_gps = gps_collection.time
    s12 = gps_collection.s12
    s18 = gps_collection.s18
    ver = [s12, s18]

    # for i_station, station in enumerate(stations):

    #     # make new plot and format 
    #     fig, ax = plt.subplots()
    #     ax.set_ylabel("<Y>: "+station + "[mm]", fontsize=14, color="green", fontweight='bold')
    #     ax.set_xlabel("GPS Time [CDT]", fontsize=14)
    #     plt.xticks(fontsize=12)
    #     plt.yticks(fontsize=12)
    #     plt.minorticks_on()
    #     ax.tick_params(axis='x', which='both', bottom=True, direction='inout')
    #     ax.tick_params(axis='y', which='both', left=True, direction='inout')
    #     plt.grid()
    #     ax.xaxis.set_major_formatter(mdate.DateFormatter('%H:%M'))
    #     fig.autofmt_xdate()
        
    #     # now, the second axes that shares the x-axis with the axs
    #     ax2 = ax.twinx()
    #     ax2.yaxis.tick_right()
    #     ax2.yaxis.set_label_position("right")
    #     ax2.set_ylabel("QHV2 [kV]", fontsize=14, color='red', fontweight='bold')
    #     ax2.xaxis.set_major_formatter(mdate.DateFormatter('%H:%M'))
        
    #     #plot data 
    #     ax.scatter(time_gps, ver[i_station], color='green', s=5.0)                          
    #     ax2.scatter(time_db, qhv, color='red', s=5.0)
         
    #     #write to disk
    #     plt.savefig("DBvsGPS_"+station+".png", dpi=300)
    #     pickle.dump(fig, open("DBvsGPS_"+station+".pickle", 'wb')) # This is for Python 3 - py2 may need `file` instead of `open`
        
    #     #plt.show()
    #     plt.cla()
    #     plt.clf()

    for i_station, station in enumerate(stations):

        # make new plot and format 
        fig, ax = plt.subplots()
        ax.set_ylabel("<Y>: "+station + " [mm]", fontsize=14, color="green", fontweight='bold')
        ax.set_xlabel("QHV [kV]", fontsize=14, color='red', fontweight='bold')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.minorticks_on()
        ax.tick_params(axis='x', which='both', bottom=True, direction='inout')
        ax.tick_params(axis='y', which='both', left=True, direction='inout')
        plt.grid()
        
        #plot data 
        ax.scatter(qhv, ver[i_station], color='green', s=12.0)                          \
        #write to disk
        plt.savefig("YvsQHV_"+station+".png", dpi=300)
        pickle.dump(fig, open("YvsQHV_"+station+".pickle", 'wb')) # This is for Python 3 - py2 may need `file` instead of `open`
    
        plt.show()



if __name__ == "__main__":
    main()
