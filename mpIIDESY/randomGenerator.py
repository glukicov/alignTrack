#
# Script to generate a plaintext file of random floating-point numbers, either 
# uniformly distributed between 0 and 1, or normally distributed with a mean 
# of 0, and a standard deviation of 1. Uses Marsenne Twister algorithm included 
# in Python's random package. User can specify seed, number of random numbers 
# to generate, precision, and output file.
#
# John Smeaton 01/02/2017
#

import sys
import getopt
import os
import random
import matplotlib.pyplot as plt


# Get system arguments, define string showing help message
argv = sys.argv[1:]
helpstring = "randomGenerator.py -s <seed> -o <output_file> -n <count> -u <uniform> -g <gaussian> -p <precision>"

# Initial values for output filename, random number seed, count of random numbers to generate, and number of decimal places to record random numbers to.
output_filename = os.path.abspath("./randoms_pre_gen.txt")
seednum = 12345678
count = 2020200
precision = 15

uniform = False
gaussian = False

# Get options, arguments
try:
    opts, args = getopt.getopt(argv, "h:s:o:n:u:g:p:", ["help", "seed", "output_file", "count", "uniform", "gaussian", "precision"])
except getopt.GetoptError:
    print helpstring
    sys.exit(2)

# Set appropriate values for seed, filename, count, and whether to generate uniform or gaussian distribution, from user arguments.
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
    elif opt in ("-u", "--uniform"):
        if arg == "True":
            uniform = True
        elif arg == "False":
            uniform = False
    elif opt in ("-g", "--gaussian"):
        if arg == "True":
            gaussian = True
        elif arg == "False":
            gaussian = False
    elif opt in ("-n", "--count"):
        count = int(arg)
    elif opt in ("-p", "--precision"):
        precision = int(arg)


if (uniform and gaussian):
    print "Please flag either uniform or gaussian, not both"
    sys.exit(2)

if ((not uniform) and (not gaussian)):
    uniform = True

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

generated_randoms = []

# Generate specified count of random numbers, writing to output file
for i in xrange(count):
    
    generated_random = 0.0;

    # Check if using uniform of gaussian distribution, writing appropriate random number to file
    if (uniform):
        generated_random = random.random()
    elif (gaussian):
        generated_random = random.normalvariate(0, 1)

    # Write to specified precision.
    output_string = "{0:." + str(precision) + "f}"
    output_file.write(output_string.format(generated_random) + "\n")
    generated_randoms.append(generated_random)

# Close output file
output_file.close()

# Sanity plot
plt.hist(generated_randoms, 100)
plt.title("Generated Randoms Sanity Histogram\n" + str(len(generated_randoms)) + " Numbers Generated, With Seed " + str(seednum))
plt.show()
