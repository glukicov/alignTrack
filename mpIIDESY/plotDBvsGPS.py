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
from scipy import stats
import numpy.polynomial.polynomial as poly
import itertools # smart lines in plotting 
# from ROOT import TH3D, TFile, TCanvas

#Pass some commands 
arg_parser = argparse.ArgumentParser(description='Input data files')
arg_parser.add_argument('--dbFile', type=str, required=True, dest='dbFile')
arg_parser.add_argument('--gpsFile', type=str, required=True, dest='gpsFile')
arg_parser.add_argument('--bins', type=int, default=np.arange(2, 100, 1), dest='bins', nargs='+')
args = arg_parser.parse_args()
dbFile = args.dbFile
gpsFile = args.gpsFile
bins = args.bins

### Constants 
db_type = [('times',np.float64),('qhv',np.float64)]
gps_type = [('times',np.float64),('s12',np.float64), ('s18',np.float64)]
stations=("S12", "S18")
# plot constants
font=14
plt.rc('xtick',labelsize=font)
plt.rc('ytick',labelsize=font)


# Data types 
db_collection_shape = ['time', 'qhv']
db_collection_reshaped_shape = ['time', 'qhv', 'eqhv']
gps_collection_shape= ['time', 's12', 's18']
gps_collection_reshaped_shape= ['time', 's12', 'es12', 's18', 'es18']
db_collection = collections.namedtuple('db_collection', db_collection_shape)
db_collection_clean = collections.namedtuple('db_collection_clean', db_collection_shape)
db_collection_reshaped = collections.namedtuple('db_collection_reshaped', db_collection_reshaped_shape)
gps_collection = collections.namedtuple('gps_collection', gps_collection_shape)
gps_collection_reshaped = collections.namedtuple('gps_collection', gps_collection_reshaped_shape)

# storage 
cors_array = [] # cors_array[i_bin][i_station]

def main():

    print("Starting plots on:", datetime.datetime.now())

    #Get data 
    print("Getting data from files...")
    db_data = getDBData(dbFile)
    gps_data = getGPSData(gpsFile)

    #Clean DB data (same start/end points)
    print("Cleaning DB data...")
    db_data_clean = cleanDB(db_data, gps_data)

    for i_bin in bins:
    
        #Rebin data 
        print("Re-binning data into",i_bin,"bins with error calculation for GPS data")
        db_data_reshaped, rescale_db = dbReshape(db_data_clean, nbins=i_bin)
        gps_data_reshaped, rescale_gps = gpsReshape(gps_data, nbins=i_bin)

        #Plot data 
        print("Plotting data..")
        cor = plotData(db_data_reshaped, gps_data_reshaped, rescale_db, rescale_gps, i_bin)
        cors_array.append(cor)
    
    plotFinal(bins, cors_array) # print cor vs bins

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
    for i_array, array in enumerate(db_collection_clean):  #lop over each array in the collection 
        array = (array[:-total_removed or None]) # remove elements 
        rescaled_array = [] 
        rescaled_array_sem = [] 
        for i_bin in range(nbins):  #loop over bins, as a range of rescaled elements 
            lower = rescale*i_bin 
            upper = rescale*(i_bin+1)
            sliced_array = array[lower:upper]
            rescaled_array.append(np.mean(sliced_array)) # take the mean of the slice (bin)
            if (i_array == 1):
                sem = stats.sem(sliced_array)
                rescaled_array_sem.append(sem)
        rescaled_data_array.append(rescaled_array)
        if (i_array == 1):
            rescaled_data_array.append(rescaled_array_sem)

    dataR = db_collection_reshaped(*rescaled_data_array)  # push arrays into the collection 
    return dataR, rescale

def gpsReshape(gps_collection, nbins=20):
    data_len = len(gps_collection[0])
    rescale = int(data_len/nbins) # elements per bin 
    total_kept = rescale * nbins # can fit in bins equally
    total_removed = data_len - total_kept # drop elements that can fit 
    rescaled_data_array =[] 
    for i_array, array in enumerate(gps_collection):  #lop over each array in the collection 
        array = (array[:-total_removed or None]) # remove elements 
        rescaled_array = [] 
        rescaled_array_sem = [] 
        for i_bin in range(nbins):  #loop over bins, as a range of rescaled elements 
            lower = rescale*i_bin 
            upper = rescale*(i_bin+1)
            sliced_array = array[lower:upper]
            rescaled_array.append(np.mean(sliced_array)) # take the mean of the slice (bin)
            # if data is the <Y> from the station 
            if (i_array==1 or i_array==2):
                sem = stats.sem(sliced_array)
                rescaled_array_sem.append(sem)
        rescaled_data_array.append(rescaled_array)
        if (i_array==1 or i_array==2):
            rescaled_data_array.append(rescaled_array_sem)

    dataR = gps_collection_reshaped(*rescaled_data_array)  # push arrays into the collection 
    return dataR, rescale
    

def plotData(db_collection, gps_collection, rescale_db, rescale_gps, i_bin):
    time_db = db_collection.time
    qhv = db_collection.qhv
    eqhv = db_collection.eqhv
    time_gps = gps_collection.time
    s12 = gps_collection.s12
    es12 = gps_collection.es12
    s18 = gps_collection.s18
    es18 = gps_collection.es18
    ver = [s12, s18]
    error = [es12, es18]

    
   # for i in range(len(qhv)):
   #      if (qhv[i] > 18.0):
   #          date = mdate.num2date(time_db[i])
   #          if (date.hour > 10):
   #              print(mdate.num2date(time_gps[i]).hour, mdate.num2date(time_gps[i]).minute, s18[i])

    #return correlation
    corr_array = []

    fig, ax = plt.subplots(2, 1, figsize=(8,10))
    
    for i_station, station in enumerate(stations):
        
        #get correlation, fit a line
        corr_matrix = np.corrcoef(qhv, ver[i_station])
        corr = np.round(corr_matrix[0][1], 3)
        x_gen = np.linspace(float(min(qhv)), float(max(qhv)), num=1000) # generate x-points for evaluation 
        coefs = poly.polyfit(qhv, ver[i_station], 1) # x1 line
        fit = poly.polyval(x_gen, coefs) # fit over generated points
        
        corr_array.append(coefs[1]) # per station 

        #plot data 
        ax[i_station].minorticks_on()
        ax[i_station].grid()
        ax[i_station].scatter(qhv, ver[i_station], color='green', label=station+":\n bins="+str(i_bin)+"\n <Y> per bin="+str(rescale_gps)+"\n QHV per bin="+str(rescale_db) )
        ax[i_station].errorbar(qhv, ver[i_station], yerr=error[i_station], xerr=eqhv, elinewidth=1, linewidth=0, capsize=2, color='green')  
        ax[i_station].plot(x_gen, fit, color="red", linestyle="--", label="Fit\n r="+str(corr)+"\n slope="+str(round(coefs[1],4))+"\n intercept="+str(round(coefs[0],4)) )
        ax[i_station].set_ylabel("<Y>: "+station + " [mm]", fontsize=font+2, fontweight='bold')
        ax[i_station].set_xlabel("QHV [kV]", fontsize=font+2, fontweight='bold')
        ax[i_station].tick_params(axis='x', which='both', bottom=True, direction='inout')
        ax[i_station].tick_params(axis='y', which='both', left=True, direction='inout')
        #get correlation, fit a line, legend 
        ax[i_station].legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 14}) # outside (R) of the plot 
        
       
    #write to disk if only passing a single bin number
    if (len(bins) < 30 ): 
        plt.tight_layout()
        plt.savefig("YvsQHV_"+str(i_bin)+".png", dpi=100)
        # pickle.dump(fig, open("YvsQHV_"+str(i_bin)+".pickle", 'wb')) # This is for Python 3 - py2 may need `file` instead of `open`

    return corr_array

def plotFinal(bins, cors_array):
   
    fig, ax = plt.subplots(2, 1, figsize=(8,10))
    
    stationM=[ [], [] ]
    for i_station, station in enumerate(stations):
 
        #plot data 
        ax[i_station].minorticks_on()
        ax[i_station].grid()
        for i_bin, the_bin in enumerate(bins):
            ax[i_station].scatter(the_bin, cors_array[i_bin][i_station], color='orange', label=station if(i_bin==0) else "")  
            stationM[i_station].append(cors_array[i_bin][i_station])
        
        mean= np.mean(stationM[i_station])
        print(station, mean)
        line =  [bins[0], mean] ,  [bins[-1], mean]
        ax[i_station].plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'black', linestyle="--")
        
        ax[i_station].set_ylabel("Slope", fontsize=font+2, fontweight='bold')
        # ax[i_station].set_ylabel(r"Correlation, $r$", fontsize=font+2, fontweight='bold'))
        ax[i_station].set_xlabel("Number of bins", fontsize=font+2, fontweight='bold')
        ax[i_station].tick_params(axis='x', which='both', bottom=True, direction='inout')
        ax[i_station].tick_params(axis='y', which='both', left=True, direction='inout')
        #get correlation, fit a line, legend 
        ax[i_station].legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 14}) # outside (R) of the plot 
        
    print("Combined mean", ( np.mean(stationM[0]) + np.mean(stationM[1]) ) /2  )
    #write to disk
    plt.tight_layout()
    plt.savefig("Final.png", dpi=100)
    #pickle.dump(fig, open("Final.pickle", 'wb')) # This is for Python 3 - py2 may need `file` instead of `open`

if __name__ == "__main__":
    main()
