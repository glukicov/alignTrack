##########################################
# Toy straw gun
# Gleb 
##########################################
import argparse, sys
import numpy as np 
import matplotlib.pyplot as plt
import subprocess, glob

#### Define constants ###
#Tracker geometry 
ModuleN=4 # TODO 8 
StrawN=8 # TODO 32
Z_seperation =  0.1 # m  / 1000 mm between gun plane and the first plane of the tracker 
X_spread = 0.015 # m / 150 mm +/- spread on the gun plane 

strawSpacing = 0.606;  # x distance between straws in a layer
layerSpacing = 0.515; # z distance between layers in a view
viewSpacing = 2.020; # z distance between views in a modules
moduleSpacing = 13.735; # z distance between modules' first layers [first layer of module 1 and first layer of module 2]
layerDisplacement = 0.303; # relative x distance between first straws in adjacent layers in a view [upstream layer is +x shifted]

hitSmearing = 0.00015 # m 150 um 

#Physics 
C = 2.99792458 * 1E8  # mm/s 
q = 1.60217 * 1E-19  # C 
B = 1.45147 # T //  assuming constant field  

P_max = 3100 # MeV

# create a defined seed RNG
np.random.seed(151214)

#Define plotting constants
COLOR="green"

### Command line input 
parser = argparse.ArgumentParser()
parser.add_argument("--nTracks", type=int, help="Number of tracks to generate")
args = parser.parse_args()

nTracks=args.nTracks

###Invoke functionality through main
def main():

    #delete files from previous run 
    clean_start()

    #Call geometry 
    geom()

    #Call gun
    x_pos, p_gen = gun(nTracks)

    #Call truth plots 
    truthPlots(x_pos, p_gen)

    #Tracking
    track()

    #Call track plots
    trackPlots()



###Define functions
def geom():
    pass
    #return stationGeom 

def gun(nTracks):
    print("Generating", nTracks, "tracks")
    
    #draw position from a UD 
    x_pos = np.random.uniform(-1, 1, nTracks)*X_spread
    
    #generate momentum (beam direction)
    p_gen = np.random.uniform(0, 1, nTracks)*P_max 

    P = [p_gen, 0] 
    
        
    return x_pos, p_gen

def truthPlots(x_pos, p_gen):
    histo(x_pos, "Generated X Position", "m")
    histo(p_gen, "Generated P", "MeV")
    pass

def track():
    pass 
    #return recoTrajectory 

def trackPlots():
    pass 

###Define helper functions
def histo(X, hist_label, units):
    # create 1 subplot
    fig0, ax0 = plt.subplots(1,1, figsize=(5,5))
    #create the histogram 
    n0,b0,p0 = ax0.hist(X, bins='auto', color=COLOR)
    ax0.set_xlabel(str(hist_label)+"/ "+str(units))
    fig0.savefig(str(hist_label.replace(" ", "_"))+".png", dpi=250)

def clean_start():
    #clean old .png files  #TODO make a function 
    glob_str="*.png"
    subprocess.Popen(["trash"] + glob.glob(glob_str))



#Define classes and data-types

#Execute code through main 
if __name__=="__main__":
    main()


