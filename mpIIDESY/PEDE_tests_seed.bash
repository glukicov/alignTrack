#!/bin/bash
clear 
echo "AlignTracker for 5000 tracks with different seeds"
#mv PEDE_Mis.txt PEDE_Mis.txt.BK
echo "Staring PEDE tests..."

x=0
#mv g_seed.txt g_seed.txt.BK
#mv u_seed.txt u_seed.txt.BK
while [ $x -le 1000 ]

do
	y=$(( ( RANDOM % 1000000 )  + 1 ))
	python randomIntGenerator.py -g True -o gaussian_ran.txt -s $y -n 80000
	echo $y >> g_seed.txt 
	sleep 2
	z=$(( ( RANDOM % 1000000 )  + 1 )) 
	python randomIntGenerator.py -u True -o uniform_ran.txt -s $z -n 5000
	echo $z >> u_seed.txt 
	sleep 1
	"./AlignTracker" d 5000
	sleep 3
	"./pede" Tracker_str.txt
	sleep 3
	python ConcatenatePEDE.py &
	sleep 0.2
	x=$(( $x + 1 ))
done

./TrackerLaunchSeed.py &
echo "PEDE Tests Complete."
echo "Producing FoM Seed plots..."

















