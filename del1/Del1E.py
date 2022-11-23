import numpy as np

from ast2000tools import utils
#utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import constants
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)

planet0_radius = system.radii[0]*1e3
planet0_mass = system.masses[0]*constants.m_sun


grav_acc = 6.67*10**(-11)*planet0_mass/planet0_radius**2
print(grav_acc*1e5*3)
print(grav_acc*1e5*3/150000)


multiplier = 6

consumption = 20*multiplier#12.07#450
init_mass = 10**5
speed_boost = 1000
thrust_force = 150000*multiplier#70974#consumption*6500


def consumed(thrust_force,consumption,init_mass,speed_boost):
      time = speed_boost*init_mass/thrust_force
      consumed = consumption*time
      time = time
      return consumed,time


r0 = planet0_radius

v_esc = np.sqrt(2*6.67*10**(-11)*planet0_mass/r0)
def gravity(r,init_mass):
    g = 6.67*10**(-11)
    return init_mass*g*planet0_mass/r**2


total_time = 0
i = 1
while total_time*speed_boost <= v_esc:
    Force = thrust_force-gravity(r0,init_mass)
    C,time = consumed(Force,consumption,init_mass,speed_boost)
    r0 = r0 + total_time*speed_boost/i
    init_mass -= C
    total_time += time
    i += 1

print(f'Remaining mass:{init_mass:.2f}kg, time:{total_time:.2f}sec')




from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)

consumption = 20*multiplier
init_mass = 10**5
speed_boost = 1000
thrust_force = 150000*multiplier

planet_x0 = np.einsum('ij->ji',system.initial_positions)[0][0]

radius =  system.radii[0]*10**3
position = np.array([utils.m_to_AU(radius)+planet_x0,0])
mission.set_launch_parameters(thrust_force,consumption,init_mass,total_time,position,0)
mission.launch_rocket(0.01)


vy0 = utils.AU_pr_yr_to_m_pr_s(np.einsum('ij->ji',system.initial_velocities)[0][1])

rot_v =  np.sqrt(2*np.pi*radius**2/(0.92282999*24*60**2))



x = utils.m_to_AU(radius+r0)+planet_x0
y = utils.m_to_AU((vy0+rot_v)*total_time)
position = np.array([x,y])
print(position)
mission.verify_launch_result(position)