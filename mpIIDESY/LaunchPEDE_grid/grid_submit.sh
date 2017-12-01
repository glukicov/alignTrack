#!/bin/bash

####Script to submit PEDE/MC iterative jobs###
# From submitting machine:
# 1) kinit -f glukicov@FNAL.GOV
# 2) ssh gm2gpvm01 [GSSAPI authentication is .ssh_config]

#--------- No action needed (e.g. already in .bash_profile)---------
# mkdir /pnfs/GM2/scratch/users/glukicov [already created]
#source /grid/fermiapp/products/common/etc/setups.sh
#setup jobsub_client
#setup fife_utils
#export SAM_EXPERIMENT=gm2
#export EXPERIMENT=gm2
#export ROLE=analysis
#---------------
# 3) jobsub_submit -N 1 -G GM2 --expected-lifetime=1h \ 
# --memory=500MB --disk=2GB \ --resource-provides=usage_model=OFFSITE \ 
# file:///gm2/app/users/glukicov/pede_test/grid_submit.sh

#Tar pede and MC dir [tracker/pede and tracker/Tracker_MC]:
# tar -zcf tracker.tar MC_pede
# scp tracker.tar 'gm2gpvm01:/gm2/app/users/glukicov/pede_test'     # TODO  how to send tar to grid?

#Data will be stored to  /pnfs/GM2/scratch/users/glukicov/ #TODO: ensure .txt files returns from the grid

#Setup 
echo "Start: "`date`
echo "Site:${GLIDEIN_ResourceName}"
echo "the worker node is " `hostname` "OS: " `uname -a`
echo "the user id is " `whoami`
echo "the output of id is " `id`

####################################################
# setup the g-2 software
####################################################
echo "Setting up g-2 software"
source /cvmfs/gm2.opensciencegrid.org/prod7/g-2/setup
setup gm2 ${RELEASE} -q prof
setup ifdh_art v1_15_05 -q e10:prof:s41 #update for v8
echo "Done .. software setup"

#################################
# set up software environments
#################################
echo "setup cvmfs common products"
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup #TODO figure out if this is needed
setup ifdhc -z /cvmfs/fermilab.opensciencegrid.org/products/common/db
setup sam_web_client
echo "done ... "

#cd $_CONDOR_SCRATCH_DIR   ? # XXX do I need this?
#echo "pwd is " `pwd`

tar -xf tracker.tar     # untar 
#Compile Fortran and C 
cd MC_pede   # C++, Fortran and Python code is here
make  # make pede 
make -f AlignTracker.mk  # make MC

#Run Fortran and C (via a python script):: stdout=open(os.devnull, 'wb') is set for MC and PEDE
python LaunchRuns.py 4 2 100000

echo "Job finished successfully on: " `date`  

#What I want to get back is a text file [line for 1 iterations with 2 runs: dm1 dm2 dm1 dm2]
# the file is (in) MC_pede/MC_pede_data.txt  - TODO job get the file

# Once I have the MC_pede_data.txt file: 
# scp 'gm2gpvm01:/gm2/app/users/glukicov/pede_test/MC_pede/MC_pede_data.txt' .
# ./plotdM.py [tested on my laptop and gm2gpvm01]
# imgcat PEDERuns.png - DONE!    

#exit 0;