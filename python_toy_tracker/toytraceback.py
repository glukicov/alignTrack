import random
import sys
import math

class Wire:
    
    # Class representing a single wire in the detector

    def __init__(self, x_rel_to_layer, x_layer_absolute, y_layer_absolute):
        
        # Finds up absolute x-position of wire, and sets other variables         
        self.x_rel_to_layer = x_rel_to_layer
        self.y_layer_absolute = y_layer_absolute
        self.absolute_x = self.x_rel_to_layer + x_layer_absolute

    def get_x_rel_to_layer(self):
        # Returns position relative to layer
        return self.x_rel_to_layer

    def get_absolute_x(self):
        # Returns absolute x-pos
        return self.absolute_x

    def get_absolute_y(self):
        # Returns absolute y-pos (same as parent layer)
        return self.y_layer_absolute

    def set_absolute_layer_x(self, x_layer_absolute):
        # Sets new position of wire, from new position of layer
        self.absolute_x = self.x_rel_to_layer + x_layer_absolute


class Layer:

    # Class representing a layer of wires

    def __init__(self, x_rel_to_plane, y_rel_to_plane, x_plane_absolute, y_plane_absolute):

        # Gets absolute position of layer
        self.x_rel_to_plane = x_rel_to_plane
        self.y_rel_to_plane = y_rel_to_plane
        self.absolute_x = self.x_rel_to_plane + x_plane_absolute
        self.absolute_y = self.y_rel_to_plane + y_plane_absolute
        
        # Creates set of wires within layer
        self.wires = [Wire(6.06 * i, self.absolute_x, self.absolute_y) for i in xrange(16)]

    def get_x_rel_to_plane(self):
        # Returns x-pos of layer relative to parent plane
        return self.x_rel_to_plane;

    def get_y_rel_to_plane(self):
        # Returns y-pos of layer relative to parent plane        
        return self.y_rel_to_plane;

    def get_wires_x(self):
        # Returns x-pos of all wires in layer
        return [wire.get_absolute_x() for wire in self.wires]
    
    def get_wires_y(self):
        # Returns y-pos of all wires in layer
        return [wire.get_absolute_y() for wire in self.wires]
    
    def get_wires(self):
        # Returns all wires in layer
        return self.wires

    def set_absolute_plane_x(self, x_plane_absolute):
        # Sets new position for layer, from new position of plane
        self.absolute_x = self.x_rel_to_plane + x_plane_absolute
        for wire in self.wires:
            wire.set_absolute_layer_x(self.absolute_x)
        

class Plane:

    # Class representing a plane of two layers of wires

    def __init__(self, x_rel_to_module, y_rel_to_module, x_module_absolute, y_module_absolute):

        # Get absolute pos of plane
        self.x_rel_to_module = x_rel_to_module
        self.y_rel_to_module = y_rel_to_module
        self.absolute_x = self.x_rel_to_module + x_module_absolute
        self.absolute_y = self.y_rel_to_module + y_module_absolute

        # Set up two layers of straws in plane
        self.layer_1 = Layer(0, 0, self.absolute_x, self.absolute_y)
        self.layer_2 = Layer(3.03, 5.15, self.absolute_x, self.absolute_y)

    def get_x_rel_to_module(self):
        # Return x-pos relative to parent module
        return self.x_rel_to_module

    def get_y_rel_to_module(self):
        # Return y-pos relative to parent module
        return self.y_rel_to_module

    def get_wires_x(self):
        # Return x-pos of all wires in plane
        return self.layer_1.get_wires_x() + self.layer_2.get_wires_x()        

    def get_wires_y(self):
        # Return y-pos of all wires in plane
        return self.layer_1.get_wires_y() + self.layer_2.get_wires_y()        

    def get_layer(self, layer_num):
        # Return layer denoted by specified number
        if layer_num == 1:
            return self.layer_1
        elif layer_num == 2:
            return self.layer_2
        
    def set_absolute_module_x(self, x_module_absolute):
        # Sets new absolute position for planes, from specified position of module
        self.absolute_x = self.x_rel_to_module + x_module_absolute
        self.layer_1.set_absolute_plane_x(self.absolute_x)
        self.layer_2.set_absolute_plane_x(self.absolute_x)


class Module:

    # Class representing a module of two planes in the detector

    def __init__(self, x_rel_to_tracker, y_rel_to_tracker, x_tracker_absolute, y_tracker_absolute):
        
        # Get absolute position of module
        self.x_rel_to_tracker = x_rel_to_tracker
        self.y_rel_to_tracker = y_rel_to_tracker
        self.absolute_x = self.x_rel_to_tracker + x_tracker_absolute 
        self.absolute_y = self.y_rel_to_tracker + y_tracker_absolute 

        # Set up two planes in module
        self.plane_1 = Plane(0, 0, self.absolute_x, self.absolute_y)
        self.plane_2 = Plane(0, 20.20, self.absolute_x, self.absolute_y)

    def get_x_rel_to_tracker(self):
        # Returns x-pos relative to parent tracker
        return self.x_rel_to_tracker;

    def get_y_rel_to_tracker(self):
        # Returns y-pos relative to parent tracker
        return self.y_rel_to_tracker;

    def get_wires_x(self):
        # Returns x-pos of all wires in module
        return self.plane_1.get_wires_x() + self.plane_2.get_wires_x()        

    def get_wires_y(self):
        # Returns y-pos of all wires in module
        return self.plane_1.get_wires_y() + self.plane_2.get_wires_y()        

    def get_plane(self, plane_num):
        # Returns plane denoted by specified number
        if plane_num == 1:
            return self.plane_1
        elif plane_num == 2:
            return self.plane_2

    def set_x_align(self, align_dist):
        # Sets x-alignment for module - gives new position to planes
        self.absolute_x = self.x_rel_to_tracker + align_dist
        self.plane_1.set_absolute_module_x(self.absolute_x)
        self.plane_2.set_absolute_module_x(self.absolute_x)
        

class Detector:    

    # Class representing the whole traceback detector

    def __init__(self):
        # Sets up two modules in the detector
        self.module_1 = Module(0, 0, 0, 0)
        self.module_2 = Module(0, 112, 0, 0)
        
    def get_wires_x(self):
        # Returns x-pos of all wires in detector
        return self.module_1.get_wires_x() + self.module_2.get_wires_x()
    
    def get_wires_y(self):
        # Returns y-pos of all wires in detector
        return self.module_1.get_wires_y() + self.module_2.get_wires_y()
        
    def get_module(self, module_num):
        # Returns module denoted by specified number
        if module_num == 1:
            return self.module_1
        elif module_num == 2:
            return self.module_2

    def set_module_x_align(self, module_num, align_dist):
        # Changes alignment of specified module.
        if module_num == 1:
            self.module_1.set_x_align(align_dist)
        elif module_num == 2:
            self.module_2.set_x_align(align_dist)


class Track:

    # Class representing a track through the detector

    def __init__(self):
        
        # Sets positions of beginning and end of track, above and below the detector.
        # Random x-positions chosen for these points
        self.y_bottom = -50
        self.y_top = 200
        self.x_bottom = random.uniform(-1, 100)
        self.x_top = random.uniform(-1, 100)

    def get_bottom_point(self):
        # Returns bottom point of track
        return [self.x_bottom, self.y_bottom]

    def get_top_point(self):
        # Returns top point of track
        return [self.x_top, self.y_top]

    def get_x_points(self):
        # Returns x-pos of beginning, end of track
        return [self.x_bottom, self.x_top]

    def get_y_points(self):
        # Returns y-pos of beginning, end of track
        return [self.y_bottom, self.y_top]


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


