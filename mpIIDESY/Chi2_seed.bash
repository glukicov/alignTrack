#!/bin/bash
clear 
echo "Bias Hunting for 50000 tracks with different seeds"
echo "Staring tests..."

x=0
while [ $x -le 2000 ]

do
	y=$(( ( RANDOM % 1000000 )  + 1 ))
	python randomIntGenerator.py -g True -o gaussian_ran.txt -s $y -n 320000
	echo $y >> g_seed.txt 
	sleep 1.4
	z=$(( ( RANDOM % 1000000 )  + 1 )) 
	python randomIntGenerator.py -u True -o uniform_ran.txt -s $z -n 60000
	echo $z >> u_seed.txt 
	sleep 0.6
	"./AlignTracker" n 20000
	sleep 8
	x=$(( $x + 1 ))
done

./TrackerLaunchSeed.py &
echo "PEDE Tests Complete."
echo "Producing FoM Seed plots..."

















