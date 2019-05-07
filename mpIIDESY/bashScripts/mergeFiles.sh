#script to merge files from a list: cat /dir/*.root > FileList.txt 
# controlled by targetFileSize [MB]
mergeDirPath=/home/gm2/glukicov/MergedInputFiles/

if [ ! -f FileList.txt ]; then
    echo "FileList.txt not found"
    return
fi

if [ -z $GM2_VERSION ]; then
    echo "GM2 software not setup"
    return 
fi

totalFiles=0 
targetFileSize=500 # MB
runningTotal=0 # per merge 
nFilesToMerge=0 # per merge 
filesToMerge="" # list passed to the gm2 process 

for file in `cat FileList.txt`; do 
  thisSize=`du -m $file | cut -f1` # in MB 
  runningTotal=$((runningTotal+thisSize)) 
  if [ -z "$filesToMerge" ]; then
    firstFile=`basename $file`
  fi
  filesToMerge=${filesToMerge}" "${file} #add to the list 
  nFilesToMerge=$((nFilesToMerge+1))   # increment 
  if [ $runningTotal -gt $targetFileSize ]; then 
    if [ $nFilesToMerge -gt 1 ]; then #check for a single large file
      echo starting a merge for $firstFile
      gm2 -c fileMerge.fcl -s $filesToMerge -o ${mergeDirPath}$firstFile
    else
      echo copying over non-merged files $filesToMerge
      cp $filesToMerge ${mergeDirPath} #for single large files simple copy to destination 
    fi
    filesToMerge="" #reset counters 
    runningTotal=0
    nFilesToMerge=0
    totalFiles=$((totalFiles+1))
  fi
done

echo done merge for 
echo totalFiles = $totalFiles

