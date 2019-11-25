# Clohessy-Wiltshire
import math
import numpy as np

z = 60											# Distance from center of station
mass = 10										# Mass of spacecraft, kg

r_earth = (np.average([6378137, 6356752]))  	# Radius of Earth, m
G = (6.67*10**-11)  							# Gravitational constant
mass_earth = (5.972*10**24)  					# Earth mass, kg
mu = G*mass_earth  								# Specific gravitational constant of Earth

at = 408000 + r_earth							# Semi-major axis of orbit (circular)
n = math.sqrt(mu/at**3)

z_pp = (n**2)*z

print(str(z_pp*mass*1000) + "mN")