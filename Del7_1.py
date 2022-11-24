#EGEN KODE
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from numba import njit
from ast2000tools import utils
from ast2000tools import constants
#utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)



### ALL DOWN TO NEXT LINE IS PART 5###
"""------------------------------------------------------------------"""

planet_positions0 = np.einsum('ij->ji',system.initial_positions)
planet_trajectories = np.load('planet_trajectories.npz')
planet_positions = np.einsum('ijk->jki',planet_trajectories['planet_positions'])
planet_velocities = np.load('velocities.npy')

print(planet_velocities.shape)
print(planet_positions.shape)


from ast2000tools import constants as con
from ast2000tools import solar_system
from ast2000tools import space_mission
from ast2000tools import shortcuts
import numpy as np
from ast2000tools import utils
#utils.check_for_newer_version()

print( utils.get_seed('Claudieg') )
seed = 93823

code_escape_trajectory =  45499
code_launch_results    =  69502



mission = space_mission.SpaceMission(seed)
system = solar_system.SolarSystem(seed)

shortcut = shortcuts.SpaceMissionShortcuts(mission, [code_escape_trajectory,code_launch_results])


thrust = 1500000 
mass_loss_rate = 200 

time =  76965*1e-4# insert the time for the spacecraft to be put on escape trajectory
height_above_suface = 10386466 # insert height above surface you want the rocket placed
direction_angle = 320# insert the angle between the x-axis and the rocket's motion
fuel_left = 42000# insert how much fuel you want for your trip

shortcut.place_spacecraft_on_escape_trajectory(
    thrust,
    mass_loss_rate,
    time,
    height_above_suface,
    direction_angle,
    fuel_left
    )



fuel_consumed, time_after_launch, pos_after_launch, vel_after_launch\
    = shortcut.get_launch_results()

mission.verify_launch_result(pos_after_launch)

mission.verify_manual_orientation(pos_after_launch,vel_after_launch,199)

travel = mission.begin_interplanetary_travel()
time, pos, vel = travel.orient()

delta_v = vel_after_launch*0.163 + np.array([0.02,0])
travel.boost(delta_v)
travel.coast(((199990)*1e-4))


delta_v = np.array([0.14360667563,-0.0499099124444])
travel.boost(delta_v)

travel.coast(((207782.1 -199990)*1e-4))
time, pos, vel = travel.orient()



def delta_v_stable_orbit(timestep,rocket_vel,rocket_pos):
    v6 = planet_velocities[timestep,6]
    p6 = planet_positions[6,timestep]
    R = rocket_pos-p6
    mu = system.masses[6]*4*np.pi**2
    #v_norm is needed velocity for stable orbit
    v_norm = np.sqrt(mu/np.linalg.norm(R))

    #v is vector orthogonal with r
    v = np.random.randn(2)
    v -= v.dot(R) * R / np.linalg.norm(R)**2
    v /= np.linalg.norm(v)
    v *= v_norm
    return v - rocket_vel + v6

delta_v = delta_v_stable_orbit(207782+76965,vel,pos)

travel.boost(delta_v)
time, pos, vel = travel.orient()

travel.coast(0.1)
time, pos, vel = travel.orient()

"""------------------------------------------------------------------"""

print('part6')


def plot_orbit(time,pos,vel):
    timestep = int(time/1e-4)
    v6 = planet_velocities[timestep,6]
    p6 = planet_positions[6,timestep]
    n = 10**3
    dt = 1e-5
    r = np.zeros((n,2))
    v = np.zeros((n,2))
    r[0] = pos - p6
    v[0] = vel - v6
    for i in range(n-1):
        R = r[i]
        a = - 4*np.pi**2*system.masses[6]*R/np.linalg.norm(R)**3
        v[i+1] = v[i] + a*dt
        r[i+1] = r[i] + v[i+1]*dt
    planet_radius = utils.m_to_AU(system.radii[6]*1e3)
    theta = np.linspace(0,2*np.pi,200)
    plt.plot(planet_radius*np.cos(theta),planet_radius*np.sin(theta),label='P6')
    plt.plot(r[:,0],r[:,1])
    plt.axis('equal')
    plt.show()


#lower velocity relative to planet by twenty percent
delta_v = -(vel-planet_velocities[int(time/1e-4),6])*0.3
travel.boost(delta_v)
travel.coast(0.1)
time, pos, vel = travel.orient()

#boost to make orbit stable
delta_v = delta_v_stable_orbit(int(time/1e-4),vel,pos)
travel.boost(delta_v)
travel.coast(0.01)
time, pos, vel = travel.orient()

#lower orbit again 
delta_v = -(vel-planet_velocities[int(time/1e-4),6])*0.3
travel.boost(delta_v)
travel.coast(0.005)
time, pos, vel = travel.orient()

#stabilize
delta_v = delta_v_stable_orbit(int(time/1e-4),vel,pos)
travel.boost(delta_v)
travel.coast(0.01)
time, pos, vel = travel.orient()

#same again
delta_v = -(vel-planet_velocities[int(time/1e-4),6])*0.3
travel.boost(delta_v)
travel.coast(0.01)
time, pos, vel = travel.orient()
delta_v = delta_v_stable_orbit(int(time/1e-4),vel,pos)
travel.boost(delta_v)
time, pos, vel = travel.orient()

#same again
delta_v = -(vel-planet_velocities[int(time/1e-4),6])*0.35
travel.boost(delta_v)
travel.coast(0.015)
time, pos, vel = travel.orient()

#repeat stable orbit a lot because it was difficult close to planet
delta_v = delta_v_stable_orbit(int(time/1e-4),vel,pos)
travel.boost(delta_v)
travel.coast(0.01)
time, pos, vel = travel.orient()


delta_v = -(vel-planet_velocities[int(time/1e-4),6])*0.15
travel.boost(delta_v)
travel.coast(0.005)
time, pos, vel = travel.orient()
delta_v = delta_v_stable_orbit(int(time/1e-4),vel,pos)
travel.boost(delta_v)
travel.coast(0.01)
time, pos, vel = travel.orient()
#last correction
delta_v = delta_v_stable_orbit(int(time/1e-4),vel,pos)
travel.boost(delta_v)
travel.coast(0.01)
time, pos, vel = travel.orient()


# travel.look_in_direction_of_planet(6)
# travel.start_video()
travel.coast(0.01)
time, pos, vel = travel.orient()
# print(pos)
# print(vel)
# travel.finish_video(filename='travel_video0.xml')

# plot_orbit(time,pos,vel)
travel.record_destination(6)
landing = mission.begin_landing_sequence()

'''-----------------------------------------'''
'''Part 6_D'''


landing.look_in_direction_of_planet(6) # Make the camera track the destination planet
time, position, velocity = landing.orient()
#landing.take_picture(filename='landing_picture.xml') 

def orbit(pos,vel):
    n = 3*10**6 #approx one orbit
    dt = 1e-2
    r = np.zeros((n,2))
    v = np.zeros(r.shape)
    r[0] = pos[0],pos[1]
    v[0] = vel[0],vel[1]
    for i in trange(n-1):
        a = -6.67*10**(-11)*system.masses[6]*constants.m_sun*r[i]/np.linalg.norm(r[i])**3
        v[i+1] = v[i] + a*dt
        r[i+1] = r[i] + v[i+1]*dt
        if (i/(7.5*10**5)).is_integer() and i != 0:
            landing.fall(7.5*10**5) 
            landing.take_picture(filename=f'LP_{int(i/(7.5*10**5))}.xml')
    plt.plot(r[:,0],r[:,1])
    theta = np.linspace(0,2*np.pi,200)
    plt.plot(system.radii[6]*1e3*np.cos(theta),system.radii[6]*1e3*np.sin(theta),label='P6')
    plt.axis('equal')
    plt.show()
#orbit(position,velocity)
'''Accordning to pictures taken. after or 7.5*10**5 if calc was wrong# seconds we are straight above a great landing position'''


def rotating_coordinates(curr_coordinates,time_elapsed):
    rot_period = system.rotational_periods[6]
    omega = 2*np.pi*rot_period * utils.day_to_s(1)

    x,y = curr_coordinates 

    R = np.linalg.norm([x,y])

    new_coordinates = np.arctan(y/x)

    new_coordinates += omega*time_elapsed

    return new_coordinates, R*np.cos(new_coordinates),R*np.sin(new_coordinates)

'''Del 7'''
'''-----------------------------------------------------'''

lander_mass = mission.lander_mass#in kg

#wait until  we are above landingsite
time,pos,vel = landing.orient()
landing.fall(7.5*10**3)
time,pos,vel = landing.orient()
landing.look_in_direction_of_planet()
landing.start_video()
#break and fall until about top of atmosphere
while np.linalg.norm(pos) -system.radii[6]*1e3 > 100000:
    landing.boost(-vel)
    landing.fall(100)
    time,pos,vel = landing.orient()

#release launcher 
landing.launch_lander(np.array([0,0,0]))

#calculate area of parachute and deploy immediately
Parachute_Area = (2*6.67*10**(-11)*lander_mass*system.masses[6]*constants.m_sun)/((system.radii[6]*1e3)**2* system.atmospheric_densities[6]*3**2)#A=2G(m_l)*M/(R^2*rho*Cd*V_safe^2)
landing.adjust_parachute_area(Parachute_Area)

from Del6_C import f2 as atmospheric_density

def gravity(r):
    return (6.67*10**(-11)*system.masses[6]*constants.m_sun*mission.lander_mass)/r**2
def thrust_force(r):
    r = np.linalg.norm(r)
    #terminal velocity
    v = np.sqrt((6.67*10**(-11)*system.masses[6]*constants.m_sun*mission.lander_mass)/(system.radii[6]*1e3)**2/Parachute_Area*atmospheric_density(r*1e-3-system.radii[6]) )
    return gravity(r) - 0.5*atmospheric_density(r*1e-3-system.radii[6])*Parachute_Area*(v**2-9)

'''ast2000 gives me slightly different numbers each time i run the code which makes it difficult to
pick a good min_height for thrust-activation
So i chose a thrustforce of 700N because the thrust_force function returns values close to 700 each time i run the code
setting a constant thrustforce which wont change each time i run the code, will make choosing the min_height a lot easier'''
F_thrust = 700#thrust_force(pos)

landing.adjust_landing_thruster(force=F_thrust,min_height=470)#(6.67*10**(-11)*system.masses[6]*constants.m_sun/(system.radii[6]*1e3)**2*(mission.lander_mass*1.155))

landing.fall(10**5)

landing.finish_video(filename='landing_video1.xml')


print(f"mass of planet: {system.masses[6]*constants.m_sun}")
print(f"radius: {system.radii[6]}")
print(f"The area is: {Parachute_Area}")
print(f"The atmospheric density is: {system.atmospheric_densities[6]}")
print(f"The period of planet 6: {system.rotational_periods[6]}")

'''
##INIT POS AND VEL PRE LANDER##
Position: (-982347, 3.714e+06, 0) m
Velocity: (159.135, -601.65, 0) m/s
V = 620m/s <- norm(velocity)

Landing thruster properties:
  Force: 700 N
  Minimum activation height: 470 m
Landing engine with thrust 700 N activated at time 814508 s.
Lander reached the surface at time 814529 s.
Successfully landed on planet 6 at time 814529 s with velocity 2.80769 m/s. Well done!
*** Achievement unlocked: Touchdown! ***
Landing site coordinates recorded:
  theta = 90 deg
  phi = 261.878 deg
'''

