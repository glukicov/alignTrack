#!/bin/bash

####Script to submit PEDE/MC iterative jobs###
# From submitting machine:
# 1) kinit -f glukicov@FNAL.GOV
# 2) ssh gm2gpvm01 [GSSAPI authentication is .ssh_config]
# 3) mkdir /pnfs/GM2/scratch/users/glukicov [already created]
# 4) source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup #set UPS area             								
# 5) setup jobsub_client 
# 6) jobsub_submit -N 2 -G GM@ --expected-lifetime=1h \ --memory=500MB --disk=2GB \ --resource-provides=usage_model=OFFSITE \ file:///gm2/app/users/glukicov/pede_test/grid_submit.sh

#Tar pede and MC dir [tracker/pede and tracker/Tracker_MC]:
# tar -zcf tracker.tar MC_pede
# scp tracker.tar 'gm2gpvm01:/gm2/app/users/glukicov/pede_test'     # TODO  how to send tar to grid?

#Data will be stored to  /pnfs/GM2/scratch/users/glukicov/ #TODO: ensure .txt files returns from the grid

#Setup 
set -x #start bash debugging at this point
echo Start `date`
echo Site:${GLIDEIN_ResourceName}
echo "the worker node is " `hostname` "OS: " `uname -a`
echo "the user id is " `whoami`
echo "the output of id is " `id`
set +x #stop bash debugging at this point

#cd $_CONDOR_SCRATCH_DIR   ? # XXX do I need this?
#echo "pwd is " `pwd`

tar xf tracker.tar     # untar 
#Compile Fortran and C 
cd MC_pede   # C++, Fortran and Python code is here
make  # make pede 
make -f AlignTracker.mk  # make MC

#Run Fortran and C (via a python script):: stdout=open(os.devnull, 'wb') is set for MC and PEDE
python LaunchRuns.py 4 2 100000

echo Job finished successfully on `date` `time` 

#What I want to get back is a text file [line for 1 iterations with 2 runs: dm1 dm2 dm1 dm2]
# the file is (in) MC_pede/MC_pede_data.txt

# Once I have the MC_pede_data.txt file: 
# scp 'gm2gpvm01:/gm2/app/users/glukicov/pede_test/MC_pede/MC_pede_data.txt' .
# ./plotdM.py [tested on my laptop and gm2gpvm01]
# imgcat PEDERuns.png - DONE!  

#exit 0;