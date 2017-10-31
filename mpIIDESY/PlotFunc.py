#!/usr/bin/env python


####################################################################
# Analytical Plotter of the residual derivatives for a range of values
#
# Created: 31 Oct 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 31 Oct 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################

import numpy as np
import matplotlib.pyplot as plt #for plotting 
from mpl_toolkits.mplot3d import Axes3D

x_range= np.arange(-40, 40, 0.1)

# dlc1 = ( c+m*z-x ) / ( np.sqrt(m*m+1) * np.absolute(c+m*z-x) ) 

# dlc2 = ( (m*m+1)*z*(c+m*z-x) - m*np.power(np.absolute(c+m*z-x),2) ) / ( np.power(m*m+1,1.5) * np.absolute(c+m*z-x)  )

# dgl1 = -1.0 * ( c+m*z-x ) / ( np.sqrt(m*m+1) * np.absolute(c+m*z-x) ) 

def graph(formula, x, name):  
    x = x  
    y = formula(x)  # <- note now we're calling the function 'formula' with x
    plt.plot(x, y)
    plt.xlabel("variable")
    plt.ylabel(name)
    plt.title(name)
    plt.show()

def dlc1(c):
    return ( c+m*z-x ) / ( np.sqrt(m*m+1) * np.absolute(c+m*z-x) )

def dlc2(c):
	return ( (m*m+1)*z*(c+m*z-x) - m*np.power(np.absolute(c+m*z-x),2) ) / ( np.power(m*m+1,1.5) * np.absolute(c+m*z-x)  )

def dgl1(c):
	return -1.0 * ( c+m*z-x ) / ( np.sqrt(m*m+1) * np.absolute(c+m*z-x) )

#Simulat. 
# m=c=z=x_range
# graph(dlc1, x_range, "dlc1: Simultaneous")
# graph(dlc2, x_range, "dlc2: Simultaneous")
# graph(dgl1, x_range, "dgl1: simultaneous")


# #Fixing slope, c, z
# m=0.0
# c=0.0
# z=31.5   
# graph(dlc1, x_range, "dlc1: Fixed m=0, c=0, z=31.5 as means")
# graph(dlc2, x_range, "dlc2: Fixed m=0, c=0, z=31.5 as means")
# graph(dgl1, x_range, "dgl1: Fixed m=0, c=0, z=31.5 as means")


# z=0.1
# m=0.1
# x=0.1   
# graph(dlc1, x_range, "dlc1(m)")
# graph(dlc2, x_range, "dlc2(m)")
# graph(dgl1, x_range, "dgl1(m)")


#3D Plots
def dlc1_3D(m, c):
    return ( c+m*z-x ) / ( np.sqrt(m*m+1) * np.absolute(c+m*z-x) )

def dlc2_3D(m ,c):
	return ( (m*m+1)*z*(c+m*z-x) - m*np.power(np.absolute(c+m*z-x),2) ) / ( np.power(m*m+1,1.5) * np.absolute(c+m*z-x)  )

def dgl1_3D(m ,c):
	return -1.0 * ( c+m*z-x ) / ( np.sqrt(m*m+1) * np.absolute(c+m*z-x) )

def graph_3D(formula, m, c, name): 
	fig = plt.figure()
	c = formula(m ,c)  # <- note now we're calling the function 'formula' with x
	ax = fig.add_subplot(111, projection='3d')
	a = [i for i in range(n) for _ in range(n)]
	b = range(n) * n
	ax.surface(a, b, c)
	ax.set_xlabel('m')
	ax.set_ylabel('c')
	ax.set_zlabel('name')
	plt.ylabel(name)
	plt.title(name)
	plt.show()

z=31.5
x=0
graph_3D(dlc1_3D, x_range, x_range, "dlc1(x=0, z=31.5)")
graph_3D(dlc2_3D, x_range, x_range, "dlc2(x=0, z=31.5)")
graph_3D(dgl1_3D, x_range, x_range, "dgl1(x=0, z=31.5)")
