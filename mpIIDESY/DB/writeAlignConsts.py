# /*
# *   Gleb Lukicov (g.lukicov@ucl.ac.uk)
# *   Created: 21 August 2019
# * /

# Read the tracker alignment offsets [mm] from a FHICL file in
# gm2tracker/align/alignmentConstants_Run1_Run2/
# and write them in the production database  

# Usage (e.g) to populate existing tables: 
# python writeAlignConsts.py --fcl=OffsetsPerModuleS12S18_Run1_Period_1_of_1_15921_17528.fcl --runStart=15921 --runEnd=19027 --iovID=1 
# Usage to create new tables: 
# python writeAlignConsts.py --createNew
# Usage to read from the DB and dump alignment constants and IoV into a file
# python writeAlignConsts.py --writeDBtoFile
# To run in a test mode with --test

import psycopg2  #db query
import json    # dbconnetion.json
import argparse # command line inputs sub
import os # MRB_SOURCE 
import re, sys # to get offsets from file
import numpy as np 
parser = argparse.ArgumentParser(description='mode')
parser.add_argument('--createNew', action='store_true', default=False,)
parser.add_argument('--fcl', type=str, default=None, help="FHICL file to load into DB")
parser.add_argument('--runStart', type=int, default=None, help="Begin of IoV")
parser.add_argument('--runEnd', type=int, default=None, help="End of IoV")
parser.add_argument('--iovID', type=int, default=None, help="The ID of this IoV")
parser.add_argument('--test', action='store_true', default=False, help="Only print out the psql commands")
parser.add_argument('--writeDBtoFile', action='store_true', default=False, help="Print the alignment tables to file")
args = parser.parse_args()
createNew=args.createNew
fcl=args.fcl
runStart=args.runStart
runEnd=args.runEnd
iovID=args.iovID
writeDBtoFile=args.writeDBtoFile
test=bool(args.test)

######Defining tracker constants ######
moduleN = 8 # 8 trackers: 
stationN = 2  # 2 stations
stations = ("12", "18")
FHICLPatchName = ["strawModuleRShift", "strawModuleHShift"] # radial and vertical offsets 
globalN = len(FHICLPatchName) # radial and vertical offsets 
FHICLServicePath = "services.Geometry.strawtracker." #full path of the tracking FHICL patch 
FHICLPath = str(os.environ['MRB_SOURCE'])+"/gm2tracker/align/alignmentConstants_Run1_Run2/" #must have local gm2tracker repo


### HELPER FUNCTIONS ######
## Read offsets from FHICL file given the FHICLPatchName
def getOffsets(FHICL_file, offset_name):
    offsets = [] #tmp storage buffer 
    for line in FHICL_file:
        if re.match(offset_name, line):
            offsets=line
            #print "offsets in func", offsets
            offsets = offsets.replace(offset_name+": [", "") 
            offsets = offsets.replace("]", "") 
            offsets = np.array([float(r) for r in offsets.split(',')])
            offsets = offsets.astype(float) 
            return offsets  
        else:
            pass

#Create new tables 
def createNewTables(cur, test):
    print "Creating new tables: 1) IoV and 2) alignment constants tables"
    
    curCommand=(
    """
    CREATE TABLE tracker_alignment_run_range_iov ( 
        iov_id  int, 
        run_start  int, 
        run_end  int 
    );
    """)
    print curCommand
    if (test==False):
        cur.execute(curCommand)
        print "tracker_alignment_run_range_iov table created"
    
    curCommand=(
    """
    CREATE TABLE tracker_alignment_offsets ( 
        iov_id  int, 
        channel  int,
        constant real 
    );
    """) 
    print curCommand
    if (test==False):
        cur.execute(curCommand)
        print "tracker_alignment_offsets table created"

def dumptoFile(cur, dumpFileName):
    f = open(dumpFileName, "w+")
    cur.execute("select * from tracker_alignment_run_range_iov;")
    rows = cur.fetchall()
    print rows
    for row in rows:
        f.write(row)
    cur.execute("select * from tracker_alignment_offsets;")
    rows = cur.fetchall()
    print rows
    for row in rows:
        f.write(row)

###========================READ FROM FHICL AND WRITE ALIGNMENT CONSTANTS===============================##
def main():
    # ###OPEN CONNECTION (as a writer!)
    dbconf = None
    with open('dbconnection_dev.json', 'r') as f: # dev 
    # with open('dbconnection_prod.json', 'r') as f: # prod 
        dbconf = json.load(f)

    # create a connection to the database using info from the json file
    cnx = psycopg2.connect(user=dbconf['user'], host=dbconf['host'], \
                           database=dbconf['dbname'], port=dbconf['port'], password=dbconf['password'])
    # query the database and obtain the result as python objects
    cur=cnx.cursor() 


    print "Running in test mode:", test
    
    #create new tables
    if(createNew):
       createNewTables(cur, test)

    #if FHICL is passed 
    if(fcl!=None):
        print "FCL locations:", FHICLPath
        # fcl_file = open(FHICLPath+"/"+fcl, "r")
        fcl_file = open(fcl, "r")
        print "Opened", fcl
        offset_input=[[ [], [] ], [ [], [] ]]

        print "runStart:", runStart, "runEnd:", runEnd, "iovID:", iovID

        #get data 
        for i_station in range(0, stationN):
            for i_global in range(0, globalN):
                offset_input[i_station][i_global] = getOffsets(fcl_file, FHICLServicePath+FHICLPatchName[i_global]+stations[i_station])
                print "S"+stations[i_station], FHICLPatchName[i_global], offset_input[i_station][i_global]
        
        response = input("Correct and write to DB? [Y=1/N=0] ")
        if(response):
            
            # write the run range 
            insCommand = "INSERT INTO tracker_alignment_run_range_iov ( iov_id, run_start, run_end ) "
            insCommand = insCommand + "VALUES ( "+ str(iovID) + "," + str(runStart) + "," + str(runEnd) + ") "
            print insCommand
            if (test==False):
                cur.execute(insCommand)
                print "IOV range recorded to DB"

            # write an offset for each of the 32 channels
            channel=1
            for i_station in range(0, stationN):
                for i_global in range(0, globalN):
                    for i_module in range(0, moduleN):

                        insCommand = "INSERT INTO tracker_alignment_offsets ( iov_id, channel, constant ) "
                        insCommand = insCommand + "VALUES ( "+ str(iovID) + "," + str(channel) + "," + str(offset_input[i_station][i_global][i_module]) + ") "
                        print insCommand
                        if (test==False):
                            cur.execute(insCommand)
                            print "IOV range recorded to DB"

                        channel+=1 

                    
        else:
            print "Exiting..."
            sys.exit()

    # read from DB 
    if(writeDBtoFile):
        dumpFileName = "tacker_offsets_iov.csv"
        dumptoFile(cur, dumpFileName)

    # close communication with the databaes before committing 
    cur.close()
    cnx.commit()
    cnx.close()

if __name__=="__main__":
    main()


# ##### Some psql shortcuts ####
# psql -U gm2_reader -h ifdbprod.fnal.gov -d gm2_online_prod -p 5452
# [connects to the Production DB, which is duplicated from the Online DB]

# \dt  - list all table
# \dn - list all schemas 
# \dt gm2tracker_sc.* - list all tables under the Tracker SC schema 
# set schema 'gm2tracker_sc';
# set schema 'gm2dq'; 

# select * from slow_control_items where name like '%HV%' limit 10;
# select * from slow_control_data where scid=946 limit 10; 

#cur.execute("select * from slow_control_data where scid=946 limit 10;")
# fetch all the rows
# rows = cur.fetchall()
# print rows


# CREATE TABLE gm2dq.tracker_hv (
#     id  SERIAL PRIMARY KEY,
#     station  smallint,
#     hv_status  BIT(64),
#     run integer,
# 	subrun integer
# );

#####Correlating Run Sunbrun with timestamp ######
#You could access the content of DAQ ODB in table "gm2daq_odb", like the following
#select json_data->'Experiment'->'Security'->'RPC hosts'->'Allowed hosts' from gm2daq_odb ;
#select run_num, json_data->'Runinfo'->'Start time' from gm2daq_odb where run_num = 8000;
# select Subrun, json_data->'Logger'->'Channels'->'0'->'Settings' from gm2daq_odb