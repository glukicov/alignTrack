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

mx_range= np.arange(-0.015, 0.015, 0.001)

def graph(formula, x1_val, name):  
    y = formula(x1_val)  # <- note now we're calling the function 'formula' with x
    plt.plot(x1_val, y)
    plt.xlabel("variable")
    plt.ylabel(name)
    plt.title(name)
    plt.show()

def dlc1(mx):
    return ( -1.0 ) / ( np.sqrt( mw*mw*my*my + mx*mx - mx*mw*my + mw*mw + 1.0 ) )

# def dlc2(m):
# 	return ( (m*m+1)*z*(c+m*z-x) - m*np.power(np.absolute(c+m*z-x),2) ) / ( np.power(m*m+1,1.5) * np.absolute(c+m*z-x)  )

# def dgl1(m):
# 	return -1.0 * ( c+m*z-x ) / ( np.sqrt(m*m+1) * np.absolute(c+m*z-x) )

#M0U0 S3 z= 5.000 x=0.182

# mw = 0.1317
# my = 0.00


# graph(dlc1, mx_range, "dlc1")

# print "dlc2(m=-0.015, c=-1.0) = ", dlc2(-0.015, -1.0)
# print "dlc2(m=0.015, c=1.0) = ", dlc2(0.015, 1.0)
# print "dlc2(m=-0.015, c=1.0) = ", dlc2(-0.015, 1.0)
# print "dlc2(m=0.015, c=-1.0) = ", dlc2(0.015, -1.0)


my=0.00034263 
mx=-0.00938614 
xt=33.8428 
yt=-0.264807

mw= -0.131652 
zw= -191.8 
xw= 35.5749 



numdlc = (my * mw * zw ) - (mx * zw) + (mw * yt) + xw - xt

print "my * mw * zw= ", my * mw * zw
print "mx * zw= ", (mx * zw)
print "(mw * yt)= ", (mw * yt)
print "xw= ", xw
print "xt= ", xt 
print "numdlc=(my * mw * zw ) - (mx * zw) + (mw * yt) + xw - xt = " , numdlc






