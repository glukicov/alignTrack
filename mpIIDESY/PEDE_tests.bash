#!/bin/bash
clear 
echo "Staring PEDE tests..."

mv PEDE_Mis.txt PEDE_Mis.txt.BK

if [ "$1" = "n" ] # first argument [n=No dca cut]
	then 
	slope=0.0001521
	intercept=0.9257
	track=( 50 250 500 1000 1500 2000 3000 4000 5000 6000 7000 8000 9000 10000 )
	timeDelay=( 1.0 1.5 1.8 2.3 3.0 3.5 4.5 6.0 7.0 9.0 10.0 11.5 12.5 14.0 )
	echo "No DCA cut"
	
else 
	slope=0.001291
	intercept=0.9730
	track=( 500 3000 8000 12000 19000 29000 40000 50000 65000 70000 80000 90000 100800 )
	timeDelay=( 1.0 1.5 2.2 2.8 3.8 5.2 6.9 8.0 11.0 11.8 13.6 14.5 16.0)
	echo "DCA cut"
fi

sleep 0.8
	
for ((i=0;i<${#track[@]};++i)); do

	echo "AlignTracker for ${track[i]}"
	# echo "Sleeping for [s]: ${timeDelay[i]} "
	"./AlignTracker" d $track[i]
	sleep ${timeDelay[i]}  #n linear complexity 
	"./pede" Tracker_str.txt
	sleep 1.5
	python ConcatenatePEDE.py &
	sleep 0.2
done

./TrackerLaunch.py &
echo "PEDE Tests Complete. "
echo "Producing FoM plots..."