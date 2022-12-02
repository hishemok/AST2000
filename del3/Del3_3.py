#EGEN KODE
import numpy as np
from tqdm import trange
from ast2000tools import utils
#utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import constants
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)


import matplotlib.pyplot as plt

def consumed(thrust_force,consumption,init_mass,speed_boost):
        time = speed_boost*init_mass/thrust_force
        consumed = consumption*time
        time = time
        return consumed,time

def gravity(r,init_mass,planet_mass):
        g = 6.67*10**(-11)
        return init_mass*g*planet_mass/r**2


''' THIS FUNCTION DID NOT WORK, BUT WE NEVER MANAGED TO FIND OUT WHY.
TH'''
def launch_from_any_planet(planet_idx,orbit_angle,angle_on_planet):

    ##PLANET POSITION AND VELOCITY IN FULL ORBIT

    n = 10**6
    v = np.zeros((n,2))
    r = np.zeros((n,2))
    dt = 6.5*1e-6

    planet_trajectories = np.load('planet_trajectories.npz')
    planet_positions = np.einsum('ijk->jki',planet_trajectories['planet_positions'])
    planet_velocities = np.load('velocities.npy')
    

    angle_to_coordinates = np.array([np.cos(orbit_angle),np.sin(orbit_angle)])
    for i in trange(len(planet_positions[0])):
        r[i] = planet_positions[0,i]
        r_unit = r[i]/np.linalg.norm(r[i])
        tol = 1e-5
        #FIND WHEN R IS APPROX AT THE WANTED COORDINATES IN RADIANS
        if abs(r_unit[0] - angle_to_coordinates[0]) < tol  and abs(r_unit[1] - angle_to_coordinates[1]) < tol:
            print(r_unit[0], r_unit[1])
            print(angle_to_coordinates[0],angle_to_coordinates[1])
            timestep = i
            print(i)
            break
    v = np.einsum('ijk->jik',planet_velocities)[planet_idx]
    r = np.einsum('ijk->jki',planet_trajectories['planet_positions'])[planet_idx]

    
    ##FIND ROCKETS POSITION ON PLANET
    planet_radius = system.radii[planet_idx]*1e3
    rocket_position = planet_radius*np.array([np.cos(angle_on_planet),np.sin(angle_on_planet)])

    ##LAUNCH FROM PLANET
    planet_mass = system.masses[planet_idx]*constants.m_sun

    #løkka her er for å se om kraften vi generer er nok for å escape planeten
    #går den noe særlig lengere enn 20 vil ikke raketten være sterk nok og vi går for fort tom for fuel
    for i in range(1,20):
        multiplier = 10+i
        m0 = 1.2*10**5

        consumption = 20*multiplier
        init_mass = m0
        speed_boost = 10
        thrust_force = 150000*multiplier

        r0 = planet_radius

        v_esc = np.sqrt(2*6.67*10**(-11)*planet_mass/r0)

        total_time = 0
        v_ = 0 
        while v_ <= v_esc:
            Force = thrust_force-gravity(r0,init_mass,planet_mass)
            C,time = consumed(Force,consumption,init_mass,speed_boost)
            v_ += Force*time/init_mass
            r0 = r0 + v_*time
            init_mass -= C
            total_time += time

        
        if total_time > 0 and init_mass > 1100:
            print(f'Remaining mass:{init_mass:.2f}kg, time:{total_time:.2f}sec')
            consumption = 20*multiplier
            init_mass = m0
            speed_boost = 10
            thrust_force = 150000*multiplier
            break


    ##CALCULATE ROCKET POSITION
    planet_x = r[timestep,0]
    planet_y = r[timestep,1]
    planet_vx = utils.AU_pr_yr_to_m_pr_s(v[timestep,0])
    planet_vy = utils.AU_pr_yr_to_m_pr_s(v[timestep,1])
    scale = r0/planet_radius
    rocket_x = utils.m_to_AU(rocket_position[0]*scale)
    rocket_y = utils.m_to_AU(rocket_position[1]*scale)

    rot_period = system.rotational_periods[planet_idx]
    rot_v = 2*np.pi*planet_radius/(rot_period*24*60**2)


    radius = planet_radius
    px = utils.m_to_AU(radius)*np.cos(angle_on_planet)+planet_x 
    py = utils.m_to_AU(radius)*np.sin(angle_on_planet)+planet_y
    position = np.array([px,py])
   
    mission.set_launch_parameters(thrust_force,consumption,init_mass,total_time,position,timestep*dt)
    mission.launch_rocket(0.01)
    total_time = 431.447 

    rot_v = 2*np.pi*radius/(0.92282999*24*60**2)

    rx = rocket_x+planet_x 
    ry = rocket_y+planet_y 
    x = rx + utils.m_to_AU(rot_v*total_time)*np.cos(np.pi/2-angle_on_planet) + utils.m_to_AU(planet_vx*total_time) 
    y = ry + utils.m_to_AU(rot_v*total_time)*np.sin(np.pi/2-angle_on_planet) + utils.m_to_AU(planet_vy*total_time)
 

    position_post_launch = np.array([x,y])
    mission.verify_launch_result(position_post_launch)

    return position_post_launch


PLANET = 0
ORBIT_ANGLE = 3.7270690724019775 #IN RADIANS
ANGLE_ON_PLANET = 3.708854754602006#IN RADIANS


rocket_position = launch_from_any_planet(planet_idx=PLANET,orbit_angle=ORBIT_ANGLE,angle_on_planet=ANGLE_ON_PLANET)







