#!/bin/bash
clear 
echo "AlignTracker for 10,000 tracks with different seeds"
#mv PEDE_Mis.txt PEDE_Mis.txt.BK
echo "Staring PEDE tests..."

x=0
#mv g_seed.txt g_seed.txt.BK
#mv u_seed.txt u_seed.txt.BK
while [ $x -le 10 ]

do
	y=$(( ( RANDOM % 1000000 )  + 1 ))
	python randomIntGenerator.py -g True -o gaussian_ran.txt -s $y -n 1618800
	echo $y >> g_seed.txt 
	sleep 3.5
	z=$(( ( RANDOM % 1000000 )  + 1 )) 
	python randomIntGenerator.py -u True -o uniform_ran.txt -s $z -n 312400
	echo $z >> u_seed.txt 
	sleep 1.5
	"./AlignTracker" n 100800
	sleep 8.0
	"./pede" Tracker_str.txt
	sleep 1.2
	python ConcatenatePEDE.py &
	sleep 0.1
	x=$(( $x + 1 ))
done

./TrackerLaunchSeed.py &
echo "PEDE Tests Complete."
echo "Producing FoM Seed plots..."

















