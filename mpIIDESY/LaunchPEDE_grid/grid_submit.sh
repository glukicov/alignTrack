#!/bin/bash

#Tar pede and MC dir [tracker/pede and tracker/Tracker_MC]:
# tar -zcf tracker.tar.gz pede/ Tracker_MC/ 

#Data will be stored to /pnfs/GM2/scratch/users/glukicov/commission/mc/<date>/

#Define constants
NRuns = 2
NIterations = 10000
NTracks = 100000 #will produce ~10,000 tracks after DCA rejection 

#Some book-keeping 
printenv
set -x #start bash debugging at this point
echo Start `date`
echo Site:${GLIDEIN_ResourceName}
echo "the worker node is " `hostname` "OS: " `uname -a`
echo "the user id is " `whoami`
echo "the output of id is " `id`
set +x #stop bash debugging at this point
cd $_CONDOR_SCRATCH_DIR
echo "pwd is " `pwd`

#Setup 
#kinit -f glukicov@FNAL.GOV ? 
source /cvmfs/gm2.opensciencegrid.org/prod7/g-2/setup
setup gm2 v7_06_06 -q prof
export PRODUCTS=$PRODUCTS:/grid/fermiapp/products/common/db
setup fife_utils

#Compile Fortran and C 
#cd tracker/pede
# make in pede/  - Makefile 
#cd tracker/Tracker_MC
# make -f AlignTracker.mk


#Run Fortran and C (via a python script)

# stdout=open(os.devnull, 'wb') is set for MC and PEDE
python LaunchRuns.py $NIterations $NRuns $NTracks


#What I want to get back is a text file [line for 1 iterations with 2 runs: dm1 dm2 dm1 dm2]
#Once I have the .txt file: plotdM.py macro will plot the results

echo Job finished successfully on `date` `time` 
exit 0;