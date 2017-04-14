#!/usr/bin/env python

import os
import string
import time
import decimal
from ROOT import TCanvas, TFile, TH1F, TTree, AddressOf, gROOT, Double
    

D=decimal.Decimal


C=[]  #stores all numbers 
F=[]


with open("C_Mp2debug_off.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers 
        # (as strings)
        number_str = line.split()
        #convert numbers to floats
        #numbers_float = [D(x) for x in numbers_str]  #map(float,numbers_str) works too
        C.append(number_str[0])
       # print numbers_float

with open("F_mp2test2_off_debug.txt") as f:
    for line in f:  #Line is a string
        #split the string on whitespace, return a list of numbers 
        # (as strings)
        number_str = line.split()
        #convert numbers to floats
        #numbers_float = [D(x) for x in numbers_str]  #map(float,numbers_str) works too
        F.append(number_str[0])
       # print numbers_float

#How much calcuations should be done by hand vs TProfile implimenation read: https://root.cern.ch/doc/master/classTProfile.html

# Make a tree
f = TFile('myTest.root','RECREATE')
t = TTree('MyTree','My test tree')

# Create a struct
gROOT.ProcessLine(\
  "struct MyStruct{\
    Int_t someInt;\
    Double_t someDouble;\
  };")

from ROOT import MyStruct

# Create branches in the tree
s = MyStruct()
t.Branch('rootInt',AddressOf(s,'someInt'),'someInt/I')
t.Branch('rootDouble',AddressOf(s,'someDouble'),'someDouble/D')

# Fill tree
for i in range(100):
  s.someInt = i
  s.someDouble = i
  t.Fill()

f.Write()
f.Close()

#c1 = TCanvas("c1","c1",900,700)
h1 = TH1F("h1", "", 10, 0, 10)
for number in C:
    print number
    h1.Fill(float(number))
#data.close()

h1.Draw("lego")



