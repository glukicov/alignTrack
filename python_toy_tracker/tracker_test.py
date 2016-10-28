import toytraceback
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
import math


def calc_approach_distance(wire, track):
    
    # Function to calculate perpendicular distance from track to wire

    # Get pos of wire
    x_0 = wire.get_absolute_x()
    y_0 = wire.get_absolute_y()
    
    # Get pos of beginning, end of track
    x_1 = track.get_bottom_point()[0]
    y_1 = track.get_bottom_point()[1]
    x_2 = track.get_top_point()[0]
    y_2 = track.get_top_point()[1]

    # Return perpendicular distance
    return abs((y_2 - y_1) * x_0 - (x_2 - x_1) * y_0 + x_2 * y_1 - y_2 * x_1) / math.sqrt((y_2 - y_1)**2 + (x_2 - x_1)**2)

    

def closest_hit_wires(detector, track):
    
    # Function to find closest approached wires in each layer in detector

    hit_wires = [] # Contains closest approached wires

    # Iterating across all layers
    for module_num in [1,2]:
        for plane_num in [1,2]:
            for layer_num in [1,2]:
                
                # Set up initial values for closest approach distance, wire in layer
                closest_approach = sys.float_info.max 
                closest_wire = None

                # Iterate across wires in layer
                for wire in detector.get_module(module_num).get_plane(plane_num).get_layer(layer_num).get_wires():
                    
                    # Calculate approach distance
                    approach = calc_approach_distance(wire, track)
                    
                    # Check if wire closer to track than previously checked ones
                    if approach < closest_approach:
                        closest_approach = approach
                        closest_wire = wire

                # Add closest approached wire to list
                hit_wires.append(closest_wire)

    # Return list of closest approached wires in each layer
    return hit_wires

# Set up detector, and get coordinates of all wires
detector = toytraceback.Detector()
x_wires = detector.get_wires_x()
y_wires = detector.get_wires_y()

# Set up track, and get pos of beginning, end points
track = toytraceback.Track()
x_track = track.get_x_points()
y_track = track.get_y_points()

# Get closest approached wires in each layer
hit_wires = closest_hit_wires(detector, track)

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
