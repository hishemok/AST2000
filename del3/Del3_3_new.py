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

def launch_from_any_planet(planet_idx,timestep,angle_on_planet):

    ##PLANET POSITION AND VELOCITY IN FULL ORBIT
    n = 10**6
    dt = 6.5*1e-6

    planet_trajectories = np.load('planet_trajectories.npz')
    planet_positions = np.einsum('ijk->jki',planet_trajectories['planet_positions'])
    planet_velocities = np.load('velocities.npy')
 
 
    v = np.einsum('ijk->jik',planet_velocities)[planet_idx]
    r = np.einsum('ijk->jki',planet_trajectories['planet_positions'])[planet_idx]
    planet_x = r[timestep,0]
    planet_y = r[timestep,1]

    planet_vx = v[timestep,0]
    planet_vy = v[timestep,1]

    # plt.plot(r[:,0],r[:,1])
    # plt.show()
    
   
    radius = system.radii[planet_idx]*1e3


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

        r0 = radius

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
            print(f'Remaining mass:{init_mass:.2f}kg, time:{total_time:.2f}sec',i)
            consumption = 20*multiplier
            init_mass = m0
            thrust_force = 150000*multiplier
            break

    
    

    rot_period = system.rotational_periods[planet_idx]
    rot_v = 2*np.pi*radius/(rot_period*utils.day_to_s(1))


    px = utils.m_to_AU(radius*np.cos(angle_on_planet)) + planet_x 
    py = utils.m_to_AU(radius*np.sin(angle_on_planet)) + planet_y
    position = np.array([px,py])
   

    mission.set_launch_parameters(thrust_force,consumption,init_mass,total_time,position,timestep*dt)
    mission.launch_rocket(0.001)
    #total_time = 431.441

#launch site position of rocket
    rx = planet_x + utils.m_to_AU(r0*np.cos(angle_on_planet))
    ry = planet_y + utils.m_to_AU(r0*np.sin(angle_on_planet))

#position of rocket after launch
    x = rx + utils.m_to_AU(rot_v*total_time)*np.cos(np.pi/2-angle_on_planet) + planet_vx*utils.s_to_yr(total_time) 
    y = ry + utils.m_to_AU(rot_v*total_time)*np.sin(np.pi/2-angle_on_planet) + planet_vy*utils.s_to_yr(total_time)


    rocket_vel_x = planet_vx + v_esc*np.cos(angle_on_planet+rot_v*total_time) 
    rocket_vel_y = planet_vx + v_esc*np.sin(angle_on_planet+rot_v*total_time)
    rocket_vel = np.array([rocket_vel_x,rocket_vel_y])


    position_post_launch = np.array([x,y])
    mission.verify_launch_result(position_post_launch)

    return position_post_launch,rocket_vel


PLANET = 0
TIMESTEP = 650
ANGLE_ON_PLANET = np.pi/2



rocket_position,rocket_velocity = launch_from_any_planet(planet_idx=PLANET,timestep=TIMESTEP,angle_on_planet=ANGLE_ON_PLANET)






