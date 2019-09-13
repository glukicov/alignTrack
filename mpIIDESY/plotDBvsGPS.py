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
# from ROOT import TH3D, TFile, TCanvas

#Pass some commands 
arg_parser = argparse.ArgumentParser(description='Input data files')
arg_parser.add_argument('--dbFile', type=str, required=True, dest='dbFile')
arg_parser.add_argument('--gpsFile', type=str, required=True, dest='gpsFile')
args = arg_parser.parse_args()
dbFile = args.dbFile
gpsFile = args.gpsFile

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
    db_data = getDBData(dbFile)
    gps_data = getGPSData(gpsFile)

    #Clean DB data (same start/end points)
    db_data_clean = cleanDB(db_data, gps_data)

    #Rebin data 
    db_data_reshaped = dbReshape(db_data_clean, nbins=2)
    gps_data_reshaped = gpsReshape(gps_data, nbins=12)

    #Plot data 
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
    # Average this every Na points so that we have nbins
    # nbins = float(nbins)
    # print("Using",nbins,"bins")
    # total_data_len = len(db_collection_clean[0])
    # print("Total of",total_data_len,"data points")
    # rescale = int(total_data_len/nbins)
    rescale = nbins
    print("rescale:",rescale)
    #Rescale data as the mean in the bin
    rescaled_data_array =[]
    for i, i_tuple in enumerate(db_collection_clean):
        print("i_tuple before: ", i_tuple[0], i_tuple[-1],"| len:", len(i_tuple)) # TODO get name
        i_tuple = np.mean(np.array(i_tuple[:int((len(i_tuple)/rescale)*rescale)]).reshape(-1,rescale), axis=1)
        print("i_tuple after: ", i_tuple[0], i_tuple[-1],"| len:", len(i_tuple)) # TODO get name
        rescaled_data_array.append(i_tuple)
    dataR = db_collection_reshaped(*rescaled_data_array)
    return dataR

def gpsReshape(gps_collection, nbins=20):
    # Average this every Na points so that we have nbins
    # nbins = float(nbins)
    # print("Using",nbins,"bins")
    # total_data_len = len(gps_collection[0])
    # print("Total of",total_data_len,"data points")
    # rescale = int(total_data_len/nbins)

    rescale = nbins
    print("rescale:",rescale)
    #Rescale data as the mean in the bin
    rescaled_data_array =[]
    for i, i_tuple in enumerate(gps_collection):
        print("i_tuple before: ", i_tuple[0], i_tuple[-1],"| len:", len(i_tuple)) # TODO get name
        i_tuple = np.mean(np.array(i_tuple[:(int(len(i_tuple)/rescale)*rescale)]).reshape(-1,rescale), axis=1)
        print("i_tuple after: ", i_tuple[0], i_tuple[-1],"| len:", len(i_tuple)) # TODO get name
        rescaled_data_array.append(i_tuple)
    dataR = gps_collection_reshaped(*rescaled_data_array)
    return dataR

def plotData(db_collection, gps_collection):
    time_db = db_collection.time
    qhv = db_collection.qhv
    time_gps = gps_collection.time[:-2 or None] # XXX super hack TODO reshape properly 
    s12 = gps_collection.s12[:-2 or None] # XXX super hack TODO reshape properly 
    s18 = gps_collection.s18[:-2 or None]   # XXX super hack TODO reshape properly 
    ver = [s12, s18]

    for i_station, station in enumerate(stations):

        # make new plot and format 
        fig, ax = plt.subplots()
        ax.set_ylabel("<Y>: "+station + "[mm]", fontsize=14, color="green", fontweight='bold')
        ax.set_xlabel("GPS Time [CDT]", fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.minorticks_on()
        ax.tick_params(axis='x', which='both', bottom=True, direction='inout')
        ax.tick_params(axis='y', which='both', left=True, direction='inout')
        plt.grid()
        ax.xaxis.set_major_formatter(mdate.DateFormatter('%H:%M'))
        fig.autofmt_xdate()
        
        # now, the second axes that shares the x-axis with the axs
        ax2 = ax.twinx()
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.set_ylabel("QHV2 [kV]", fontsize=14, color='red', fontweight='bold')
        ax2.xaxis.set_major_formatter(mdate.DateFormatter('%H:%M'))
        
        #plot data 
        ax.scatter(time_gps, ver[i_station], color='green', s=5.0)                          
        ax2.scatter(time_db, qhv, color='red', s=5.0)
         
        #write to disk
        plt.savefig("DBvsGPS_"+station+".png", dpi=300)
        pickle.dump(fig, open("DBvsGPS_"+station+".pickle", 'wb')) # This is for Python 3 - py2 may need `file` instead of `open`
        
        #plt.show()
        plt.cla()
        plt.clf()

    for i_station, station in enumerate(stations):

        # make new plot and format 
        fig, ax = plt.subplots()
        ax.set_ylabel("<Y>: "+station + "[mm]", fontsize=14, color="green", fontweight='bold')
        ax.set_xlabel("QHV [kV]", fontsize=14, color='red', fontweight='bold')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.minorticks_on()
        ax.tick_params(axis='x', which='both', bottom=True, direction='inout')
        ax.tick_params(axis='y', which='both', left=True, direction='inout')
        plt.grid()
        
        #plot data 
        ax.scatter(qhv, ver[i_station], color='green', s=5.0)                          \
        #write to disk
        plt.savefig("YvsQHV_"+station+".png", dpi=300)
        pickle.dump(fig, open("YvsQHV_"+station+".pickle", 'wb')) # This is for Python 3 - py2 may need `file` instead of `open`
        
        plt.show()
        plt.cla()
        plt.clf()



    # plot_array = [] # keep arrays in scope
    # tfile = TFile("DBvsGPS.root", "new")
    # for i_station, station in enumerate(stations):
    #     print(station)
    #     plot = TH3D(str(station), str(station), 100, 737230.0, 737230.9, 24, 9.0, 23.0, 14, -7.0, 9.0)
    #     # gps and db have different shapes 
    #     for i_entry in range(len(time_db)):
    #         plot.Fill(time_db[i_entry], qhv[i_entry], 0)
    #     for i_entry in range(len(time_gps)):
    #         plot.Fill(time_gps[i_entry], 0, ver[i_station][i_entry])
    #     plot_array.append(plot)
    #     plot.Draw()
    #     plot.Write()
    # tfile.Close() 


if __name__ == "__main__":
    main()
