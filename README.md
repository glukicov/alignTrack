# alignTrack

Welcome to this collaborative repo for UCL g-2 team!  

The file mp2tst.bin [default binary file] is produced by mptest2.f90 [lines 98, 354- 381]. 

Here are the instructions to get the code on working on DAQ1: 

0) scl enable devtoolset-3 python27 bash

1)  git clone https://github.com/glukicov/alignTrack.git
to get the latest code from our repository 

2) cd alignTrack/mpIIDESY

3) make
to build the pede executable 

4) test that it works by pede -t 
(should give a terminal output [last 2 lines]:
 Millepede II-P ending   ... Mon Dec 12 12:31:15 2016 
 Peak dynamic memory allocation:    0.100512 GB

5) make -f handler.mk
This build my handler code from MilleHandler.cpp (requires C++11) 
It has root capabilities for future testing/integration 

6) Run the above with  ./MilleHandler
This produces test.bin and test.root

7) To check the binary file do python readMilleBinary.py "FILENAME" "N of line"
e.g. python readMilleBinary.py test.bin -1 [reads our binary for all lines] 
e.g. python readMilleBinary.py mp2tst.bin 2 [reads default binary for 2 lines] 
