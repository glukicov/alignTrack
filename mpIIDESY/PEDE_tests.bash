#!/bin/bash
clear 
echo "Staring PEDE tests..."

mv PEDE_Mis.txt PEDE_Mis.txt.BK

for i in 50 250 500 1000 1500 2000 3000 4000 5000 6000 7000 8000 9000 10000  # no cut 
#for i in 500 3000 8000 12000 19000 29000 40000 50000 65000 70000 85000 90000 101000 # DCA cut

do
	sleep 0.2
	echo "AlignTracker for $i"
	"./AlignTracker" d $i
	sleep 25
	"./pede" Tracker_str.txt
	sleep 1.5
	python ConcatenatePEDE.py &
	sleep 0.2
done
./TrackerLaunch.py &
echo "PEDE Tests Complete."
echo "Producing FoM plots..."

















