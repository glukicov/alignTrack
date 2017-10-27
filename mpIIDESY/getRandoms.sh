#!/bin/sh


rand=$(expr $1 \* 16)
norm=$(expr $1 \* 4)

python randomIntGenerator.py -g True -o gaussian_ran.txt -s 987654321 -n $rand

python randomIntGenerator.py -u True -o uniform_ran.txt -s 123456789 -n $norm