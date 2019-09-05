#!/bin/env python

"""

  ConDBJsonFileWriter.py
	Script designed to convert fhicl files to json files.

  author: T. Walton (twalton@fnal.gov)

  July 2019

  updated by G. Lukicov (g.lukicov@ucl.ac.uk) 4 Sep 2019
  --added tracker alignment functionality 

"""


import os, sys, math, subprocess, re
import numpy as np

from optparse import OptionParser

"""
helper classes
"""
_parser = OptionParser()



"""
script option parser
"""
def JsonParserHelpOpts() :
    _parser.add_option("--fcl", dest="fcl", default=None, type="string", help="The input FHiCL file to convert to json file. [Required]")
    _parser.add_option("--file", dest="file", default="filename", type="string", help="The based name of the output json file name.")
    _parser.add_option("--path", dest="path", default="outfilePath", type="string", help="The path of the output json files.")
    _parser.add_option("--scratch", dest="scratch", default=False, action="store_true", help="Write to /pnfs scratch area and overwrite path directory.")
    _parser.add_option("--iov", dest="iov", default=None, type="string", help="The interval of validity. [Required]")
    _parser.add_option("--ana", dest="ana", default=None, type="string", help="The type of analysis. [Required]")    


"""
make string
"""
def MakeString(lst) :
    lstring = ""
    for l in lst :
        lstring += "\"%s\"," % l
    return lstring[:-1]

"""
get the filename 
"""
def GetFilename(opts,identifier) :
    jsonFilename = opts["file"].replace(".json","")
    topDir       = opts["path"]

    if opts["scratch"] :
      subdir = ""
      iovs   = opts["iov"].split(":")
      for iov in iovs :
          subdir += iov.capitalize()

      user   = str(os.environ.get('USER'))
      topDir = "/pnfs/GM2/scratch/users/%s/JSONFiles/%s/ID%s" % (user,subdir,identifier)

    if not os.path.exists(topDir) : 
       os.makedirs(topDir)

    jsonFilename = "%s/%s%s.json" % (topDir,jsonFilename,identifier)
    return jsonFilename 


"""
create a json per IOV
"""
def CreateJsonFile(opts,run) :
    
    # chose identider based on iov type 
    if opts["iov"] == "run" :
      identifier   = int(run)
    if opts["iov"] == "run_range" :
      identifier   = str(run[0])+"_"+str(run[1])
      
    jsonFilename = GetFilename(opts,identifier) 
    jsonFile     = open(jsonFilename,'w')
    jsonFile.write("{\n")
       
    #write to the json file based on the iov type 
    if opts["iov"] == "run" :
      jsonFile.write("\"interval_of_validity\" : { \"run\" : %s },\n\n" % run)

    # if run_range is used, run is an array: run = [run_start, run_end]
    if opts["iov"] == "run_range" :
       jsonFile.write( " \"interval_of_validity\" : { \"run_start\" : "+str(run[0])+", \"run_end\" : "+str(run[1])+" },\n\n " )

    return jsonFile

"""
create json file header
"""
def CreateJsonFileHeader(jsonFile,folderName,columnString,typeString) :
    jsonFile.write("\"%s\" : \n" % folderName)
    jsonFile.write("{\n")
    jsonFile.write("\t\"columns\" : [ %s ],\n" % columnString)
    jsonFile.write("\t\"types\" : [ %s ],\n" % typeString)


"""
write the json file constants
"""
def WriteJsonFileConstants(jsonFile,values,writeComma=True) :
    if not jsonFile.closed :

       for idx, value in enumerate(values) :
           svalue      = [str(i) for i in value]
           deliminator = "," if idx != len(values)-1 else " ]\n"
           if idx == 0 :
              jsonFile.write("\t\"values\" : [ [%s]%s\n" % (','.join(svalue),deliminator))   
           else :
              jsonFile.write("\t\t[%s]%s\n" % (','.join(svalue),deliminator)) 

       if writeComma :
          jsonFile.write("},\n\n")
       else :
          jsonFile.write("}\n\n")


"""
write the json file constants as an array
Input: 2D array 
"""
def WriteJsonFileConstants2DArray(jsonFile,values,writeComma=True) :
    if not jsonFile.closed :
       jsonFile.write(" \t\"values\" : [ ")
       #write a dimension of the array, with elements comma separated (repr)
       inner, outer, elements = np.shape(values) # return array shape 
       total_length = inner + outer
       i_length=0 # checking if reached the last dimension 
       for i_inner in range(inner):
          for i_outer in range(outer):
              i_length+=1 
              jsonFile.write( str(repr(values[i_inner][i_outer])))
              if (i_length != total_length):
                jsonFile.write("\n") # write new line for all but last array

       # write the final closing bracket 
       jsonFile.write(" ]\n")
       
       if writeComma :
          jsonFile.write("},\n\n")
       else :
          jsonFile.write("}\n\n")
 
"""
closing the json file
"""
def CloseJsonFile(opts,jsonFile,folderStatusName,run) :
    if not jsonFile.closed :
       cstring = ""
       tstring = ""
       values  = None

       if opts["iov"] == "run" :
          cstring = "\"run\", \"status\""
          tstring = "\"int\", \"int\""
          values  = [ [run,1] ]

       # if run_range is used, run is an array: run = [run_start, run_end]
       # as passed in the CreateJsonFiles() function 
       if opts["iov"] == "run_range" :
          cstring = " \"run_start\", \"run_end\",  \"status\" "
          tstring = " \"int\", \"int\", \"int\" "
          values  = [ [run[0], run[1], 1] ]

       CreateJsonFileHeader(jsonFile,folderStatusName,cstring,tstring)
       WriteJsonFileConstants(jsonFile,values,False) 

       jsonFile.write("}\n")
       jsonFile.close()

"""
Tracker: read offsets from a FHICL file given the FHICLPatchName
"""
def getOffsets(FHICL_file):
    
    # return data structure offsets[station][globalDoF]
    # globalDoF : radial or vertical shift 
    offsets_array = [ [ [], [] ], [ [], [] ] ]
    
    #Defining tracker constants ##
    stations = ("12", "18")
    FHICLPatchName = ["strawModuleRShift", "strawModuleHShift"] # radial and vertical offsets 
    FHICLServicePath = "services.Geometry.strawtracker." #full path of the tracking FHICL patch 

    # loop over stations and globalDoF 
    # to match the line are sore offsets in the array 
    for line in FHICL_file:
      for i_station in range(0, len(stations)):
        for i_global in range(0, len(FHICLPatchName)):
          offset_name = FHICLServicePath+FHICLPatchName[i_global]+stations[i_station]
          if re.match(offset_name, line):
              # for the matched line strip all but the real numbers 
              offsets=line
              offsets = offsets.replace(offset_name+": [", "") 
              offsets = offsets.replace("]", "") 
              offsets = [float(r) for r in offsets.split(',')]
              offsets_array[i_station][i_global]=offsets

    FHICL_file.close()
    return offsets_array 

"""
primary function to create json files
"""
def CreateJsonFiles(opts) :
    fclFile  = open(opts["fcl"],'r') 
    
    if opts["ana"] == "kicker" :
       folderStatusName = "kicker_dqc_status"
       foldername = "kicker_dqc"
       columns    = [ "kickerAmplitudeUpperBound", "kickerAmplitudeLowerBound", "kickerTimingUpperBound", "kickerTimingLowerBound", "kickerFWHMUpperBound", "kickerFWHMLowerBound", "zeroCrossingCut" ]
       types      = [ "real", "real", "real", "real", "real", "real", "int" ]

       cstring    = MakeString(columns)
       tstring    = MakeString(types)

       run        = None
       values     = None 
       jsonFile   = None

       # loop over file lines
       for fclLine in fclFile.readlines() :
           line   = fclLine.strip()

           # create a new file 
           if line.find("run") != -1 and line.find(":") != -1 :
              values   = [ [0.,0.,0.,0.,0.,0.,0.], [0.,0.,0.,0.,0.,0.,0.], [0.,0.,0.,0.,0.,0.,0.] ] 
              run      = line[:line.find(":")].replace("run","").strip() 
              jsonFile = CreateJsonFile(opts,run)
              CreateJsonFileHeader(jsonFile,foldername,cstring,tstring)
           elif line.find("}") != -1 :
              if jsonFile != None :
                 WriteJsonFileConstants(jsonFile,values)
                 CloseJsonFile(opts,jsonFile,folderStatusName,run)
           else :
             key  = -1
             if line.find("Kicker1") != -1 :
                key = 0
             elif line.find("Kicker2") != -1 :
                key = 1
             elif line.find("Kicker3") != -1 :
                key = 2
             else :
                continue

             item = -1 
             if line[:line.find("Kicker")] in columns :
                item = columns.index(line[:line.find("Kicker")]) 
             else :
                continue

             if types[item] == "real" :
                values[key][item] = float(line[line.find(":")+1:])
             elif types[item] == "int" :
                values[key][item] = int(line[line.find(":")+1:])

    # if alignment FHICL file is passed write a JSON file with constants
    # as an array of 8 offsets; 4 channels: S12_rad, S12_ver, S18_rad, S18_rad
    # total of 32 offsets per IoV (run_range)
    elif opts["ana"] == "tracker_align" :
    
       folderStatusName = "tracker_align_status"
       foldername = "tracker_align"
       columns    = [ "offsets" ]
       types      = [ "real[8]" ] 

       cstring    = MakeString(columns)
       tstring    = MakeString(types)

       run_range  = [ None, None ] # run_start and run_end 
       values     = None 
       jsonFile   = None

       # this returns a 2x2 array: values[station][DoF]
       values = getOffsets(fclFile)
       # re-open the file and get the run_start and run_end 
       fclFile  = open(opts["fcl"],'r') 
       lines=fclFile.readlines()
       run_range[0] = [int(s) for s in lines[0].split() if s.isdigit()][0]  # run_start
       run_range[1] = [int(s) for s in lines[1].split() if s.isdigit()][0]  # run_end
       
       # write the header once 
       jsonFile = CreateJsonFile(opts,run_range)
       CreateJsonFileHeader(jsonFile,foldername,cstring,tstring)
       # 2D array is being loped over its dimensions inside the function  
       WriteJsonFileConstants2DArray(jsonFile,values)
       # write the tail and close the file 
       CloseJsonFile(opts,jsonFile,folderStatusName,run_range)
    ###end of tracker_align statement 

    else :
      fclFile.close()
      sys.exit( "The analysis is unsupported. Please update the script to account for your FHiCL file format." )


   


"""
main function
"""
def main() :

    # Get the options
    if sys.argv[1] == "--help" or sys.argv[1] == "-h" :
       JsonParserHelpOpts()
       _parser.parse_args( "--help".split() )

    else :
       JsonParserHelpOpts()
       (opts, args) = _parser.parse_args()
       opts = vars(opts)

       if opts["fcl"] == None :
          sys.exit("The input fhicl file is required.")

       if opts["iov"] == None :
          sys.exit("The interval of validity must be specify.")

       if opts["ana"] == None :
          sys.exit("The analysis type must be specify.")     

    # Create Json files
    CreateJsonFiles(opts)
     


if __name__ == "__main__" :
   main()

