#!/usr/bin/env python

#TODOs in C/F: "back" calculations 


import os
#import numpy as np
#import matplotlib.pyplot as plt
import string
import time
import decimal
D=decimal.Decimal



C=[]
F=[]


os.system("make -f MakeCTest.mk")
time.sleep(1)
print "Hit Enter is compilation was successful"
test = raw_input()
if (test != ""): 
    sys.exit()

os.system("gfortran -o F_test F_test.f90")
time.sleep(1)
print "Hit Enter is compilation was successful"
if (test != ""): 
    sys.exit()


os.system("./C_test")
time.sleep(2)

os.system("./F_test")
time.sleep(2)



with open("C_P_test.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers 
        # (as strings)
        numbers_str = line.split()
        #convert numbers to floats
        numbers_float = [float(x) for x in numbers_str]  #map(float,numbers_str) works too
        C.append(numbers_float)
       # print numbers_float

with open("F_P_test.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers 
        # (as strings)
        numbers_str = line.split()
        #convert numbers to floats
        numbers_float = [float(x) for x in numbers_str]  #map(float,numbers_str) works too
        F.append(numbers_float)
       # print numbers_float


print "     floats          doubles             floats(MUL)            doubles(MUL)"
print "C: ", C
print "F: ", F
print ""

for i in range(0,4):

    str_C= str(D(C[i][0]))
    str_F= str(D(F[i][0]))
    
    list_C=[]
    list_F=[]

    for digit in str_C:
        if (digit == "-"):
            continue
        if (digit == "."):
            continue
        list_C.append(digit)
     

    for digit in str_F:
        if (digit == "-"): 
            continue
        if (digit == "."):
            continue
        list_F.append(digit)

    for n in range(0, len(list_C)-1):
        if (list_C[n] != list_F[n]):
            
            if (i==0): 
                print "The difference in floats C - F is at the", n+1,"th digit"
                break

            if (i==1): 
                print "The difference in doubles C - F is at the", n+1,"th digit"
                break

            if (i==2): 
                print "The difference in floats C - F is with MUL=1000 is at the", n+1,"th digit"
                break

            if (i==3): 
                print "The difference in doubles C - F with MUL=1000 is at the", n+1,"th digit"
                break


