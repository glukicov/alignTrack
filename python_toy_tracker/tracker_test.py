import toytraceback
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
import math


# Set up detector, and change alignment of first module
detector = toytraceback.Detector()
detector.set_module_x_align(1, 10.0)

# Get wire coordinates
x_wires = detector.get_wires_x()
y_wires = detector.get_wires_y()

# Set up track, and get pos of beginning, end points
track = toytraceback.Track()
x_track = track.get_x_points()
y_track = track.get_y_points()

# Get closest approached wires in each layer
hit_wires = toytraceback.closest_hit_wires(detector, track)

# Get pos of all closest approached wires
x_hits = [wire.get_absolute_x() for wire in hit_wires]
y_hits = [wire.get_absolute_y() for wire in hit_wires]

# Plot points where wires located, where wires hit, and line showing track
plt.scatter(x_wires, y_wires)
plt.scatter(x_hits, y_hits, color='green')
plt.plot(x_track, y_track, 'r')

# Label axes, and show plot
plt.axis('equal')
plt.xlabel("x-position / mm")
plt.ylabel("y-position / mm")
plt.show()
