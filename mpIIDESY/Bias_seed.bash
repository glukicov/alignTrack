#!/bin/bash
clear 
echo "Bias Hunting for 500 tracks with different seeds"
echo "Staring tests..."

x=0
while [ $x -le 8 ]

do
	y=$(( ( RANDOM % 1000000 )  + 1 ))
	python randomIntGenerator.py -g True -o gaussian_ran.txt -s $y -n 800000
	echo $y >> g_seed.txt 
	sleep 1.8
	z=$(( ( RANDOM % 1000000 )  + 1 )) 
	python randomIntGenerator.py -u True -o uniform_ran.txt -s $z -n 150000
	echo $z >> u_seed.txt 
	sleep 0.8
	"./AlignTracker" n 100000
	sleep 18
	x=$(( $x + 1 ))
done

./TrackerLaunchSeed.py &
echo "PEDE Tests Complete."
echo "Producing FoM Seed plots..."

















