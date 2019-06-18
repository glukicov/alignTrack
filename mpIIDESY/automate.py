import subprocess
import numpy as np
import glob

x_min = np.arange(300, 3000, 300)
x_max = np.arange(600, 3300, 300)
y_min = np.array([-25, -8, -2,   -2,  -2, -15, -15, -15, -65 ])
y_max =  np.array([5,   3,   3,   3,   3,   8,    8,  8,   10 ])


for i_step in range(0, len(x_min)):
    subprocess.call([ "python3", "/Users/gleb/software/alignTrack/mpIIDESY/CurveProfile.py", "--mode=graph", "--coord=radial", "--y_min="+str(y_min[i_step]), "--y_max="+str(y_max[i_step]),"--x_min="+str(x_min[i_step]),"--x_max="+str(x_max[i_step]) ])
    subprocess.call(["mv", "radial_graph_Pos_vs_mom_timeCut.png", "radial_graph_Pos_vs_mom_timeCut_"+str(i_step)+".png"])
