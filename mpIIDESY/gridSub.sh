# Sean's grid file 
#!/bin/bash                                                                                                                                                 
​
WORKINGDIR=$PWD
​
FCLDIR=/gm2/app/users/sbfoster/TrackingFolder/studies/trackerblock/gridOfTargets/tracking
fclFile=RunTrackingSim.fcl
OUTDIR=/pnfs/GM2/scratch/users/sbfoster/Projects/TrackerAlignmentBlock/studies/gridOfTargets/tracking
PRODDIR=/gm2/app/users/sbfoster/TrackingFolder/gm2Dev_v9_29_00/srcs/gm2analyses/ProductionScripts/produce/
RINGSIMDIR=/pnfs/GM2/scratch/users/sbfoster/Projects/TrackerAlignmentBlock/studies/gridOfTargets/ringsim
​
#make a file list in ring sim directory (make sure there is only one directory file)
if [ $(find $RINGSIMDIR/* -maxdepth 0 -type d | wc -l) == 1 ]; then
    RINGSIMDIR2=$(find $RINGSIMDIR/* -maxdepth 0 -type d)
    #delete existing file list if it exists
    echo $RINGSIMDIR2
    rm -rf $RINGSIMDIR2/FileList.txt
    #make new file list
    ls $RINGSIMDIR2/data/*.root > $RINGSIMDIR2/FileList.txt
else
    echo "Found zero or more than one directories in the ringsim output directory."
    exit 1
fi
​
cd $PRODDIR
​
./gridSetupAndSubmitGM2Data.sh --tag=tracker --mc --full --fhiclFile=$FCLDIR/$fclFile --gasGun --localArea --output-dir=$OUTDIR --njobs=100 --lifetime=8h --schemas=None --offsite=no --noifdh_art --input-filelist=$RINGSIMDIR2/FileList.txt
​
cd $WORKINGDIR