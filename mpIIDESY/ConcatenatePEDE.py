#!/usr/bin/python

####################################################################
# Concatenates PEDE results for comparison
#
# 
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################

# Getting constants from MC
with open("Tracker_p_constants.txt") as f:
     for line in f:  #Line is a string
        number_str = line.split()
        moduleN=int(number_str[0])
        trackN=int(number_str[4])
        

print "Parameters from Simulation:"
print "moduleN= ",moduleN
print "trackN= ",trackN


label = []
mis = []
error =[]

# Reading 1 set of PEDE results
with open("millepede.res") as f:
    first_line = f.readline() # skip header
    for line in f:  #Line is a string
        number_str = line.split()
        label.append(int(number_str[0]))
        mis.append(float(number_str[1]))
        error.append(float(number_str[4]))

f = open('PEDE_Mis.txt', 'a')
for i in range (0, moduleN):
    f.write(str(label[i]))
    f.write(" ")
    f.write(str(mis[i]))
    f.write(" ")
    f.write(str(error[i]))
    f.write(" ")
f.write(str(trackN))
f.write("\n")
f.close()  

print "Misalignments, errors, labels and track # appended to PEDE_Mis.txt"