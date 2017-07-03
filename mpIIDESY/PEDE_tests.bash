#!/bin/bash
clear 
echo "Staring PEDE tests..."

mv PEDE_Mis.txt PEDE_Mis.txt.BK

for i in 50 250 500 1000 1500 2000 3000 4000 5000 6000 7000 8000 9000 10000

do
	sleep 0.2
	echo "AlignTracker for $i"
	"./AlignTracker" d $i
	sleep 3.5
	"./pede" Tracker_str.txt
	sleep 2
	python ConcatenatePEDE.py &
	sleep 0.2
done
./TrackerLaunch.py &
echo "PEDE Tests Complete."
echo "Producing FoM plots..."

















