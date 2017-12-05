#!/bin/sh


rand=$(expr $1 \* 16)
norm=$(expr $1 \* 4)

seed_uni=$2
# seed_uni=$(expr $2 \*5)
# seed_gaus=$(expr $2 \*3)
seed_gaus=$2 

python randomIntGenerator.py -g True -o gaussian_ran.txt -s $seed_gaus -n $rand

python randomIntGenerator.py -u True -o uniform_ran.txt -s $seed_uni -n $norm