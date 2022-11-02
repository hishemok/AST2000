#EGEN KODE

import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from numba import njit
from ast2000tools import utils
from ast2000tools import constants
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)



planet_positions0 = np.einsum('ij->ji',system.initial_positions)

planet_trajectories = np.load('planet_trajectories.npz')
planet_positions = np.einsum('ijk->jki',planet_trajectories['planet_positions'])
planet_velocities = np.load('velocities.npy')

print(planet_velocities.shape)
print(planet_positions.shape)
'''
distance close enough to destination planet
l = norm(r)*sqrt(mp/(10ms))
'''


###1year
'''
one orbit is equal to 8.6143yrs
p0 = planet_positions[0]
for i in range(1,p0.shape[0]):
    if p0[i,1] > 0 and p0[i-1,1] < 0:
        print(i,i*1e-4)
        break
'''


from ast2000tools import constants as con
from ast2000tools import solar_system
from ast2000tools import space_mission
from ast2000tools import shortcuts
import numpy as np
from ast2000tools import utils
utils.check_for_newer_version()

print( utils.get_seed('Claudieg') )
seed = 93823# insert student's seed number
# sjekk at dette er hva dere får ^

code_escape_trajectory =  45499
code_launch_results    =  69502

"""------------------------------------------------------------------"""

mission = space_mission.SpaceMission(seed)
system = solar_system.SolarSystem(seed)

shortcut = shortcuts.SpaceMissionShortcuts(mission, [code_escape_trajectory,code_launch_results])

################################################################
#            PLACE SPACECRAFT ON ESCAPE TRAJECTORY             #
################################################################
#                   |      For Part 3      |
#                   ------------------------

# OBS!!! Here the angle is in degrees!
# Copy and paste from here, but choose only one of the two parts (A or
# B) before copying the rest. Parts A and B depends on what the
# student(s) actually need this shortcut for. The parts are:

# SHORTCUT A:
# use this if the student(s) only need shortcut for generalised launch,
# i.e., were able to launch their spacecraft in Part 1.

# SHORTCUT B:
# use this if the student(s) needed Shortcut B in the two previous
# shortcuts.
"""
------------------------------------------------------------------------
"""




"""
Documentation
place_spacecraft_on_escape_trajectory(
    rocket_thrust, rocket_mass_loss_rate, time height_above_surface,
    direction_angle, remaining_fuel_mass):

------------------------------------------------------------------------
place_spacecraft_on_escape_trajectory() places the spacecraft on an
escape trajectory pointing directly away from the home planet.

Parameters
----------
rocket_thrust  :  float
    The total thrust of the rocket, in NEWTONS.

rocket_mass_loss_rate  :  float
    The total mass loss rate of the rocket, in KILOGRAMS PER SECOND.

time  :  float
    The time at which the spacecraft should be placed on the escape
    trajectory, in YEARS from the initial system time.

height_above_surface  :  float
    The heigh above the home planet surface to place the spacecraft, in
    METERS (after launch).

direction_angle  :  float
    The angle of the direction of motion of the spacecraft with respect
    to the x-axis, in DEGREES.

remaining_fuel_mass  :  float
    The mass of fuel carried by the spacecraft after placing it on the
    escape trajectory, in KILOGRAMS.

Raises
------
RuntimeError
    When none of the provided codes are valid for unlocking this method.
------------------------------------------------------------------------

"""

# ----------
# Shortcut A
# ----------

thrust = 1500000 # insert the thrust force of your spacecraft here
mass_loss_rate = 200 # insert the mass loss rate of your spacecraft here

# ----------
# Shortcut B (if the student(s) need compute_engine_performance())
# ----------


"""
GROUP TEACHER OF AST2000:
------------------------------------------------------------------------
"""
# Use shortcut B of the shortcut compute_engine_performance for thrust
# and mass_loss_rate. The students that need this part should already
# have this when they used the shortcut for Part 1.B. 
"""
------------------------------------------------------------------------
"""


# choose these values freely, but they should be relevant to where you
# want to go, e.g., if you want to travel outwards of your solar system,
# don't let the direction angle be 0 if you are launching from
# coordinates close to (-x, 0), as this will send you in the opposite
# direction), and vice versa if your destination is a planet closer to
# your sun
'''Timestep from hohmann transfer function "Del5_1"
Angle launches about same direction as planet trajectory at timestep
Height above surface is approx where we are after escape velocity is reached'''
time =  76965*1e-4# insert the time for the spacecraft to be put on escape trajectory
height_above_suface = 10386466 # insert height above surface you want the rocket placed
direction_angle = 320# insert the angle between the x-axis and the rocket's motion
fuel_left = 30000# insert how much fuel you want for your trip

shortcut.place_spacecraft_on_escape_trajectory(
    thrust,
    mass_loss_rate,
    time,
    height_above_suface,
    direction_angle,
    fuel_left
    )

"""
GROUP TEACHER OF AST2000:
------------------------------------------------------------------------
"""
# In order to get the position of the spacecraft after launch, the
# students will have to use this shortcut from Part 1 of the project,
# i.e., this should not affect their score for Part 3.
"""
------------------------------------------------------------------------
"""

fuel_consumed, time_after_launch, pos_after_launch, vel_after_launch\
    = shortcut.get_launch_results()
# fuel_consumed is None because place_spacecraft_on_escape_trajectory()
# don't actually launch the rocket, but just place the rocket in correct
# position and with the correct velocity.
mission.verify_launch_result(pos_after_launch)

mission.verify_manual_orientation(pos_after_launch,vel_after_launch,199)

travel = mission.begin_interplanetary_travel()
time, pos, vel = travel.orient()

delta_v = vel_after_launch*0.163 + np.array([0.02,0])
travel.boost(delta_v)
#time, pos, vel = travel.orient()

    
travel.coast(((199990)*1e-4))

#time, pos, vel = travel.orient()

delta_v = np.array([0.14360667563,-0.0499099124444])#np.array([0.14,0])
travel.boost(delta_v)
#time, pos, vel = travel.orient()

travel.coast(((207782.1 -199990)*1e-4))

time, pos, vel = travel.orient()
#We are now within the distance l and are going to try to stay in a stable orbit 



'''
tried finding how far i was from my destination planet at end of travel
used that information to boost the rocket in the right direction. Thats why the prev boost is so specific.

print('mypos = ',pos)
print('planetpos', planet_positions[6,207782+76965])
print('delta_r',planet_positions[6,207782+76965]-pos)
print('delta_r_norm',utils.AU_to_m(np.linalg.norm(planet_positions[6,207782+76965]-pos)))
'''

def delta_v_stable_orbit(timestep,rocket_pos,rocket_vel):
    p6 = planet_positions[6,timestep]
    v6 = planet_velocities[timestep,6]
    rocket = rocket_pos
    vel = rocket_vel
    print('vel,v6',vel,v6)
    orbit_vel_norm = np.sqrt(4*np.pi**2*system.masses[6]/np.linalg.norm(rocket))
    pos_vec = p6-rocket

    #create random vector before making it orthogonal on pos_vec
    vel_orbit = np.random.randn(2) 
    vel_orbit -= vel_orbit.dot(pos_vec) * pos_vec / np.linalg.norm(pos_vec)**2

    #make the orbit velocity vector the right "size", 
    # so norm(vel_orbit) = orbit_vel_norm which is the needed veloctiy to orbit the planet
    vel_orbit *= orbit_vel_norm/np.linalg.norm(vel_orbit)

    #return the boost needed to leave the rocket in a stable orbit around planet 6
    delta_v = vel_orbit - vel + v6
    print('vel_orb',vel_orbit)
    return delta_v



# travel.look_in_direction_of_planet(6)
# travel.start_video()
delta_v = delta_v_stable_orbit(207782+76965,pos,vel)

travel.boost(delta_v)
time, pos, vel = travel.orient()
print('p6',planet_positions[6,207782+76965])
print('v6',planet_velocities[207782+76965,6])

travel.coast(.01)
time, pos, vel = travel.orient()
print(planet_positions[6,207782+76965 + 10**2])

# travel.finish_video(filename='travel_video.xml')
