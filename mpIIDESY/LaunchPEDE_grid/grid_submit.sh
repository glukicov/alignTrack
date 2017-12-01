#!/bin/bash

#Tar pede and MC dir [tracker/pede and tracker/Tracker_MC]:
# tar -zcf tracker.tar MC_pede
# scp tracker.tar 'gm2gpvm01:/gm2/app/users/glukicov/pede_test'
# scp grid_submit.sh 'gm2gpvm01:/gm2/app/users/glukicov/pede_test'

#Data will be stored to /pnfs/GM2/scratch/users/glukicov/commission/mc/<date>/

#Setup 
#kinit -f glukicov@FNAL.GOV ? 
source /cvmfs/gm2.opensciencegrid.org/prod7/g-2/setup
setup gm2 v7_06_06 -q prof
export PRODUCTS=$PRODUCTS:/grid/fermiapp/products/common/db
setup fife_utils

tar xf tracker.tar
#Compile Fortran and C 
cd MC_pede
make
make -f AlignTracker.mk

#Run Fortran and C (via a python script)
# stdout=open(os.devnull, 'wb') is set for MC and PEDE
python LaunchRuns.py 4 2 100000

echo Job finished successfully on `date` `time` 
#What I want to get back is a text file [line for 1 iterations with 2 runs: dm1 dm2 dm1 dm2]
#Once I have the .txt file: plotdM.py macro will plot the results

# the file is (in) MC_pede/MC_pede_data.txt
# scp 'gm2gpvm01:/gm2/app/users/glukicov/pede_test/MC_pede/MC_pede_data.txt' .
# ./plotdM.py
# imgcat PEDERuns.png - DONE!  

#exit 0;