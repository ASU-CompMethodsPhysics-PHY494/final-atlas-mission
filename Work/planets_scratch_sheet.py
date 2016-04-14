import numpy as np

# mass of solar system bodies in kg
body_mass = {'Sun': 1.989e30, 'Venus': 4.87e24, 'Earth': 5.97e24, 
             'Mars': 0.642e24, 'Jupiter': 1898e24, 'Saturn': 568e24, 
             'Titan': 0.1345e24}

# distance from sun in km
distance = {'Venus': 108.2e6, 'Earth': 149.6e6, 'Mars': 227.9e6,
             'Jupiter': 778.6e6, 'Saturn': 1433.5e6}

initial_radial_position = {'Venus': 14.22, 'Earth': 219.09, 'Mars': 229.93, 
                           'Jupiter': 172.26, 'Saturn': 252.77}

def polar_to_cartesian(initial_radial_position, distance):
    """converts polar coordinates to cartesian coordinates
    Input:
    ------
    initial_radial_position
        initial angle in degrees
    distance
        distance from sun to planet (used for conversion calculation)
        
    Output:
    -------
    x, y
        x and y coordinates where sun is at origin
    """
    rad = deg*(np.pi/180)
    x = distance*np.cos(rad)
    y = distance*np.sin(rad)
    return x, y
    