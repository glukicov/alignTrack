############
# A single script that runs residual plotter, PEDE, and alignment plotter
# 
# Assumed to be run from a folder with :
# 1) TrackerAlignemnt.root (alignment ana file from RunAlignmentPlots.fcl)
# 2) PEDE data file from Alignment producer (e.g. trackerAlignment_S12.fcl)
#
# Created: 28 March 2019 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 28 March 2019 by Gleb
##############
import subprocess

print("Removing old plots")
subprocess.call(["trash", "PEDE_Results_S12.png"])
subprocess.call(["trash", "PEDE_Results_S18.png"])
subprocess.call(["trash", "FoM_Res_S12.png"])
subprocess.call(["trash", "FoM_Res_S18.png"])

print("Calling residual plotter")
subprocess.call(["python3", "/Users/gleb/software/alignTrack/mpIIDESY/GetRes.py"])

print("Calling pede in the local dir")
subprocess.call(["/Users/gleb/software/alignTrack/PEDE/pede", "SteeringFile.txt"]) 

print("Calling the alignement plotting script")
subprocess.call(["python3", "/Users/gleb/software/alignTrack/mpIIDESY/RobustTrackerLaunch.py"])

print("Opening new plots")
subprocess.call(["/Users/gleb/.iterm2/imgcat", "FoM_Res_S12.png"])
subprocess.call(["/Users/gleb/.iterm2/imgcat", "FoM_Res_S18.png"])
subprocess.call(["/Users/gleb/.iterm2/imgcat", "PEDE_Results_S12.png"])
subprocess.call(["/Users/gleb/.iterm2/imgcat", "PEDE_Results_S18.png"])