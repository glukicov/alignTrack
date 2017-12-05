#!/usr/bin/env python


####################################################################
# Analytical Plotter of the residual derivatives for a range of values
#
# Created: 31 Oct 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 31 Oct 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################

import numpy as np
import matplotlib.pyplot as plt #for plotting 
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

m_range= np.arange(-0.015, 0.015, 0.001)
c_range= np.arange(-1.0, 1.0, 0.001)



def graph(formula, x1_val, name):  
    y = formula(x1_val)  # <- note now we're calling the function 'formula' with x
    plt.plot(x1_val, y)
    plt.xlabel("variable")
    plt.ylabel(name)
    plt.title(name)
    plt.show()

def dlc1(m):
    return ( c+m*z-x ) / ( np.sqrt(m*m+1) * np.absolute(c+m*z-x) )

def dlc2(m):
	return ( (m*m+1)*z*(c+m*z-x) - m*np.power(np.absolute(c+m*z-x),2) ) / ( np.power(m*m+1,1.5) * np.absolute(c+m*z-x)  )

def dgl1(m):
	return -1.0 * ( c+m*z-x ) / ( np.sqrt(m*m+1) * np.absolute(c+m*z-x) )

#M0U0 S3 z= 5.000 x=0.182
z=5.5149999
x=-0.1209999
c=1.0
graph(dlc1, m_range, "dlc1")
graph(dlc2, m_range, "dlc2")

print "dlc2(m=-0.015, c=-1.0) = ", dlc2(-0.015, -1.0)
print "dlc2(m=0.015, c=1.0) = ", dlc2(0.015, 1.0)
print "dlc2(m=-0.015, c=1.0) = ", dlc2(-0.015, 1.0)
print "dlc2(m=0.015, c=-1.0) = ", dlc2(0.015, -1.0)


