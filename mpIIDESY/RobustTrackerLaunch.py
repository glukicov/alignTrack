####################################################################
# Sanity plots for Tracker Alignment.  
# FoM Plots for comparison of actual misalignment vs PEDE results 
# Used as a final step of FixedFoM.py script 
# Can deal with removed modules from tacking
#
# Created: 29 November 2018 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 29 November 2018 by Gleb
#####################################################################

import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import argparse, sys


#df = pd.read_csv("/Users/gleb/software/gm2code/PY/Python-Machine-Learning-Second-Edition/Chapter02/iris.data", header=None) #offline
df = pd.read_csv("", header=None) #offline

print df
