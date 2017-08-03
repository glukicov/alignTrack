#!/bin/bash
clear 
echo "Bias Hunting for 5000 tracks with different seeds"
echo "Staring tests..."

x=0
while [ $x -le 500 ]

do
	y=$(( ( RANDOM % 1000000 )  + 1 ))
	python randomIntGenerator.py -g True -o gaussian_ran.txt -s $y -n 80000
	echo $y >> g_seed.txt 
	sleep 0.5
	z=$(( ( RANDOM % 1000000 )  + 1 )) 
	python randomIntGenerator.py -u True -o uniform_ran.txt -s $z -n 15000
	echo $z >> u_seed.txt 
	sleep 0.5
	"./AlignTracker" d 5000
	sleep 3.0
	x=$(( $x + 1 ))
done

./TrackerLaunchSeed.py &
echo "PEDE Tests Complete."
echo "Producing FoM Seed plots..."

















