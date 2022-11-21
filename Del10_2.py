#EGEN KODE
import numpy as np
from ast2000tools import utils
from ast2000tools import constants
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)

sigma = constants.sigma
c = constants.c
k_B = constants.k_B
G = constants.G
m_p = constants.m_p
# H_m = 1.6735575*1e-27
# HE_m = 6.6464731*1e-27
# mu = (H_m*0.75+HE_m*0.25)/H_m

#def sun info
R_sun = constants.R_sun
L_sun = constants.L_sun
M_sun = constants.m_sun

#our star
star_radius = system.star_radius*1e3
star_temperature = system.star_temperature
star_mass = system.star_mass


#A
rho = star_mass*M_sun/(4/3*np.pi*star_radius**3)
#use mu = 1
T_c = star_temperature + 2*np.pi/3 * G*rho*m_p*star_radius**2/k_B
print(T_c)
#15987854.884627542
#16mill Kelvin

