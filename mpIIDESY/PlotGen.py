#!/usr/bin/python

# TODOs:
# 1. Implement where python takes in # straw, layer, modules, views, and reads and stores sensibly TODO 

####################################################################
# Sanity plots for AlginTracker
#
# 
#
# Created: 16 May 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 16 May 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################

import numpy as np 
import matplotlib.pyplot as plt #for plotting 
import itertools
import csv

#as quadruplets x0,y0,x1,y1 for Generated tracks 
x0=[] 
z0=[] 
x1=[] 
z1=[] 

#as quintuplets Ixs0, Ixs1, Ixs2, Ixs3, Izs for straw Ideal geometry 
Ixs0=[]
Ixs1=[]
Ixs2=[]
Ixs3=[] 
Izs=[]


#as quintuplets Mxs0, Mxs1, Mxs2, Mxs3, Mzs for straw Misaligned geometry 
Mxs0=[]
Mxs1=[]
Mxs2=[]
Mxs3=[] 
Mzs=[]

x_fit=[]

w, h = 4, 4;
Mis = [[0 for x in range(w)] for y in range(h)] 
Ideal = [[0 for x in range(w)] for y in range(h)] 

#TODO need a better reading method: read by columns
#Read file and store in lists 
i=0
with open("Tracker_d_geom.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
		Ideal[i][0]=(float(number_str[0]))
		Ideal[i][1]=(float(number_str[1]))
	  	Ideal[i][2]=(float(number_str[2]))
	  	Ideal[i][3]=(float(number_str[3]))
	  	Izs.append(float(number_str[4]))
	  	i=i+1



#Read file and store in lists 
with open("Tracker_p_gen.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number_str = line.split()    
        x0.append(float(number_str[0]))
        z0.append(float(number_str[1]))
        x1.append(float(number_str[2]))
        z1.append(float(number_str[3]))

with open("Tracker_p_fit.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number_str = line.split()    
        x_fit.append(float(number_str[0]))
        
 #Read file and store in lists
i=0
with open("Tracker_d_mis.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number_str = line.split()    
        Mis[i][0]=(float(number_str[0]))
        Mis[i][1]=(float(number_str[1]))
        Mis[i][2]=(float(number_str[2]))
        Mis[i][3]=(float(number_str[3]))
        Mzs.append(float(number_str[4]))
        i=i+1

# print Mis[0]
# print Mis[1]
# print Mis[2]
# print Mis[3]
# print Mzs

plt.figure(1)
plt.subplot(221)
for i in range(0, len(x0)): 
	dataM = [[x0[i],z0[i]], [x1[i],z1[i]]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(dataM, 2))),
	    color = 'brown', marker = '*')

for i in range(0, 4):
	plt.plot(Mis[0][i], Mzs[0], color="green", marker = "o")
	plt.plot(Mis[1][i], Mzs[1], color="green", marker = "o")
	plt.plot(Mis[2][i], Mzs[2], color="green", marker = "o")
	plt.plot(Mis[3][i], Mzs[3], color="green", marker = "o")


axes = plt.gca()
axes.set_xlim([-1,3])
axes.set_ylim([-3,28])
plt.ylabel("z [cm]")
plt.xlabel("x [cm]")
plt.title("Mis. (Real) Geom")


plt.subplot(222)
for i in range(0, len(x_fit)): 
	dataI = [[x_fit[i], z0[i]], [x_fit[i], z1[i]]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(dataI, 2))),
	    color = 'purple', marker = 'x')


for i in range(0, 4):
	plt.plot(Ideal[0][i], Izs[0], color="yellow", marker = "o")
	plt.plot(Ideal[1][i], Izs[1], color="yellow", marker = "o")
	plt.plot(Ideal[2][i], Izs[2], color="yellow", marker = "o")
	plt.plot(Ideal[3][i], Izs[3], color="yellow", marker = "o")

axes = plt.gca()
axes.set_xlim([-1,3])
axes.set_ylim([-3,28])
plt.ylabel("z [cm]")
plt.xlabel("x [cm]")
plt.title("Ideal Geom.")

plt.show()
