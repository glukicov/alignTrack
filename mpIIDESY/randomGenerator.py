#
# Script to generate a plaintext file of random floating-point numbers between 0
# and 1. Uses Marsenne Twister algorithm included in Python's random package. 
# User can specify seed, number of random numbers to generate, and output file.
#
# John Smeaton 27/01/2017
#

import sys
import getopt
import os
import random


# Get system arguments, define string showing help message
argv = sys.argv[1:]
helpstring = "randomGenerator.py -s <seed> -o <output_file> -n <count>"

# Initial values for output filename, random number seed, and count of random numbers to generate
output_filename = os.path.abspath("./randoms_pre_gen.txt")
seednum = 12345678
count = 2020200

# Get options, arguments
try:
    opts, args = getopt.getopt(argv, "h:s:o:n:", ["help", "seed", "output_file", "count"])
except getopt.GetoptError:
    print helpstring
    sys.exit(2)

# Set appropriate values for seed, filename, and count from user arguments.
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print helpstring
        sys.exit(2)
    elif opt in ("-s", "--seed"):
        try:
            seednum = int(arg)
        except Exception:
            print "Invalid seed - must be castable to integer"
            sys.exit(2)
    elif opt in ("-o", "--output_file"):
        output_filename = os.path.abspath(arg)

# Create directory for output file, if it doesn't already exist
if not os.path.exists(os.path.dirname(output_filename)):
    try:
        os.makedirs(os.path.dirname(output_filename))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise


# Seed random number generator
random.seed(seednum)

# Open output file
output_file = open(output_filename, "w")

# Generate specified count of random numbers, writing to output file
for i in xrange(count):
    output_file.write(str(random.random()) + "\n")

# Close output file
output_file.close()
