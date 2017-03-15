import numpy as np
import matplotlib.pyplot as plt

filename = "uniform_read.txt"

read_values = []

with open(filename) as histfile:
    for line in histfile.readlines():

        items = line.split()

        read_values.append(float(items[0]))
        

plt.hist(read_values, 50)
plt.show()
