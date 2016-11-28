#!/bin/sh
scl enable python27 "
source /home/smeaton/python27-hep/bin/activate
python /home/smeaton/g-2/python_toy_tracker/multi_fit.py > /home/smeaton/g-2/python_toy_tracker/log.txt"
