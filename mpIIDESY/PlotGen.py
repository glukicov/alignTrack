#!/usr/bin/python

# TODOs:
# 1. 

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

arr = []
arr.append([])

#TODO need a better reading method: read by columns
#Read file and store in lists 
i=0
with open("Tracker_d_geom.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
		arr[i].append(float(number_str[0]))
		arr[i].append(float(number_str[1]))
	  	arr[i].append(float(number_str[2]))
	  	arr[i].append(float(number_str[3]))
	  	arr[i].append(float(number_str[4]))
	  	i=i+1
    
        
   

	

#Read file and store in lists 
with open("Tracker_d_geom.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number_str = line.split()    
        x0.append(float(number_str[0]))
        z0.append(float(number_str[1]))
        x1.append(float(number_str[2]))
        z1.append(float(number_str[3]))

 #Read file and store in lists 
with open("Tracker_d_mis.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number_str = line.split()    
        Mxs0.append(float(number_str[0]))
        Mxs1.append(float(number_str[1]))
        Mxs2.append(float(number_str[2]))
        Mxs3.append(float(number_str[3]))
        Mzs.append(float(number_str[4]))

fig=plt.figure()
ax=fig.add_subplot(111)
for i in range(0, len(x0)): 
	data = [[x0[i],z0[i]], [x1[i],z1[i]]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(data, 2))),
	    color = 'brown', marker = '*')
for i in range(9, len(Ixs0)):
	plt.plot(Ixs0[i], Izs[0], color="green", marker = "o")
	plt.plot(Ixs0[i], Izs[0], color="green", marker = "o")
	plt.plot(Ixs0[i], Izs[0], color="green", marker = "o")
	plt.plot(Ixs0[i], Izs[0], color="green", marker = "o")

axes = plt.gca()
axes.set_xlim([-1,3])
axes.set_ylim([-3,28])
plt.title("Generated Tracks")

#plt.show()
