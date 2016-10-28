#
# Simple script to show array of wires in two modules, with track through them.
# May be adapted for graphical output of fitting program.
#

import matplotlib.pyplot as plt
import numpy as np
import random

# To contain x, y coordinates of wires
wire_points = []

# Displacements between layers, modules
layer_disp = [0, 20.20]
module_disp = [0, 112]


for module in xrange(2):

    for layer in xrange(2):
        
        for wire in xrange(7):

            # Get y-position of wire based on layer and module, then x-position based on wire number
            y = layer_disp[layer] + module_disp[module]
            x = wire * 3.03
        
            # Adjust y-position of wire based on layer and wire number, to produce staggered arrangement 
            if (layer == 1) & (wire % 2 == 0):
                y += 5.15               
            elif (layer == 0) & (wire % 2 != 0):
                y += 5.15
                
            # Add this wire's coordinates to array
            wire_points.append([x,y])


# Get arrays of x-, y-positions for wires, in order to plot
x_wires = np.array(wire_points)[:,0]
y_wires = np.array(wire_points)[:,1]

# y-coordinates for origin, ends of track
top_y = 200
bottom_y = -50

# Generate random x-coordinates for ends of track
top_x = random.uniform(-1, 20)
bottom_x = random.uniform(-1, 20)

# Get array of x-, y-positions for ends of track, for plotting
x_track = [bottom_x, top_x]
y_track = [bottom_y, top_y]

# Plot points where wires located, and line showing track
plt.scatter(x_wires, y_wires)
plt.plot(x_track, y_track, 'r')

# Label axes, and show plot
plt.axis('equal')
plt.xlabel("x-position / mm")
plt.ylabel("y-position / mm")
plt.show()
