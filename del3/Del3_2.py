import numpy as np
import matplotlib as pyplot
from ast2000tools import constants
from ast2000tools import utils

seed = utils.get_seed('Claudieg')#Claudieg

from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)

sigma = constants.sigma
T = system.star_temperature


rs = system.star_radius*10**3 #radius star
R = system.radii[0]*10**3 #radius planets

A_p = 2*np.pi*R**2 #area of planet facing sun
A_s = 4*np.pi*rs**2 #area of star

L = 4*np.pi*rs**2*sigma*T**4 #luminosity L


for planet_idx in range(system.number_of_planets):
    r = utils.AU_to_m(np.linalg.norm(np.einsum('ij->ji',system.initial_positions)[planet_idx]))
    F = (rs/r)**2*sigma*T**4         #calculates flux received per unit area by planet
    print(f"The flux received by planet {planet_idx} per unit area is: {F}")

    area = 40/(F*0.12)    #calculates the area of the solar panel, 40W, 12 percent efficiency
    print(f"the minimum area of the solar panel to power the spacecraft at planet {planet_idx} must be {area}m^2")
