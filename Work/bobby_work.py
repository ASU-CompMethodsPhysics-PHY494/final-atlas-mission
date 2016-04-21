import numpy as np

class body:
	"""defines a new class for large bodies in our system
	Input:
	------
	mass
		mass of the body measured in kg
	radius
		linear distance from the Sun (center of system) measured in km
	angle	
		angular position measured in degrees
	v0
		tangential velocity measured in km/day
	"""
	def __init__(self,mass,radius,angle,v0):
		self.mass = mass
		self.radius = radius
		self.angle = angle
		self.v0 = v0

Sun = body(1.989e30,0,0,0)
Venus = body(4.87e24,108.2e6,14.22,1.9414e7)
Earth = body(5.97e24,149.6e6,219.09,3.15533e7)
Mars = body(0.642e24,227.9e6,229.93,5.93568e7)
Jupiter = body(1898e24,778.6e6,172.26,7.05542e7)
Saturn = body(568e24,1433.5e6,252.77,1.308528e8)
Titan = body(0.1345e24,1434.7e6,252.77,4.81248e5)


def polar_to_cartesian(initial_deg, distance):
    """converts polar coordinates to cartesian coordinates
    Input:
    ------
    initial_deg
        initial angle in degrees
    distance
        distance from sun to planet (used for conversion calculation)
        
    Output:
    -------
    x, y
        x and y coordinates where sun is at origin
    """
    rad = initial_deg*(np.pi/180)
    x = distance*np.cos(rad)
    y = distance*np.sin(rad)
    return np.array([x, y])

initial_positions = polar_to_cartesian(initial_deg['Venus'], distance['Venus'])
print(initial_positions)
