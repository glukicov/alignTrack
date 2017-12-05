#!/bin/bash

####Script to submit PEDE/MC iterative jobs###
# 1) kinit -f glukicov@FNAL.GOV
# 2) ssh gm2gpvm01 [GSSAPI authentication is in .ssh_config]

#Tar pede and MC dir:
# tar -zcf tracker.tar MC_pede
#Data will be stored to  /pnfs/GM2/scratch/users/glukicov/ #TODO: return .txt file from the grid via ifdh cp

#Setup 
echo "Start: "`date`
echo "Site:${GLIDEIN_ResourceName}"
echo "the worker node is " `hostname` "OS: " `uname -a`
echo "the user id is " `whoami`
echo "the output of id is " `id`

#################################
# set up software environments
#################################
echo "setup cvmfs common products"
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup #TODO figure out if this is needed
setup ifdhc -z /cvmfs/fermilab.opensciencegrid.org/products/common/db
setup sam_web_client

echo "done setting up cvmfs common products "

####################################################
# setup the g-2 software
####################################################
#echo "Setting up g-2 software"
source /cvmfs/gm2.opensciencegrid.org/prod7/g-2/setup
setup gm2 v7_06_06 -q prof
#source /cvmfs/gm2.opensciencegrid.org/prod/g-2/setup
#source /cvmfs/oasis.opensciencegrid.org/gm2/prod7/g-2/setup
#source /cvmfs/gm2.opensciencegrid.org/prod7/g-2/setup
#setup gm2 v7_06_06 -q prof
#setup ifdh_art v1_15_05 -q e10:prof:s41 #update for v8
#echo "Done setting up g-2 software setup"

ifdh cp /pnfs/GM2/scratch/users/glukicov/tracker.tar ./tracker.tar
echo "ls:" `ls`
tar -xvzf tracker.tar     # untar 

#Compile Fortran and C 
cd MC_pede   # C++, Fortran and Python code is here
make  # make pede 
make -f AlignTracker.mk  # make MC

#Run Fortran and C (via a python script):: stdout=open(os.devnull, 'wb') is set for MC and PEDE
# Launch for #iterations #runs #tracks #seed TODO 
python LaunchRuns.py 1 2 100000

ifdh cp MC_pede_data.txt /pnfs/GM2/scratch/users/glukicov/pede_results/MC_one_pede_data_${PROCESS}.txt

echo "Job finished successfully on: " `date`  

#What I want to get back is a text file [line for 1 iterations with 2 runs: dm1 dm2 dm1 dm2]
# the file is (in) MC_pede/MC_pede_data.txt

# Once I have the MC_pede_data.txt file: 
# scp 'gm2gpvm01:/gm2/app/users/glukicov/pede_test/MC_pede/MC_pede_data.txt' .
# ./plotdM.py [tested on my laptop and gm2gpvm01]
# imgcat PEDERuns.png - DONE!    

#Get job log: jobsub_fetchlog -G gm2 --jobid=
exit 0