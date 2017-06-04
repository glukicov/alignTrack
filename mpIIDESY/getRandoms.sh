#!/bin/sh



python randomIntGenerator.py -g True -o gaussian_ran.txt -s 987654321 -n $1

python randomIntGenerator.py -u True -o uniform_ran.txt -s 123456789 -n $1