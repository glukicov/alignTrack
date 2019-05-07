iterDirPath=/home/gm2/glukicov/ana_aM6 # now this can be incremented for automating iterations 
fhiclName=RunTrackingPlots_mdc1.fcl # constant  

if [ -f STOP ]; then
  echo "STOP file found. Not running"
  return
fi

if [ -z $GM2_VERSION ]; then
    echo "GM2 software not setup"
    return 
fi

if [ $# -ne 1 ]; then
  echo "Exactly one argument must be supplied. $# were found."
  return
fi

# Get filename from argument 1
filename=$1
echo "Running process on file ${filename}..."

# Chop after gm2tracker...
fileid=${filename##*gm2tracker_particle_gun_full_}
# Chop before .root
fileid=${fileid%%.root*}

# Output file name
outputFile=${iterDirPath}/gm2tracker_ana_${filename##*gm2tracker_particle_gun_full_}

# Check that output file doesn't exist already
if [ -f $outputFile ]; then
  echo "Output file (${outputFile}) already exists.  Not running..."
  return
fi 

# Make a directory that we'll run this job in
mkdir $fileid
cd $fileid

# Run this job in its own directory
gm2 -c ${iterDirPath}/${fhiclName} -s $filename -T $outputFile

cd ../