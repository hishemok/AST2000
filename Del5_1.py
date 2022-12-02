#EGEN KODE

import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from ast2000tools import utils
from ast2000tools import constants
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)

planet_positions_at_t0 = np.einsum('ij->ji',system.initial_positions)

planet_trajectories = np.load('planet_trajectories.npz')
planet_positions = np.einsum('ijk->kji',planet_trajectories['planet_positions'])
planet_velocities = np.load('velocities.npy')


''' this function is not the one we used to pull the most accurate data from.
This one is more of a test function and a starting simulation to see what we would have to take into consideration.
The only thing we really used in this file is the hohmann transer angle at the bottom
The function we chose to use to simulate landing is in Del5_1_2f.py'''
def rocket_trajectory_2(position,velocity,time):
    dt = 1e-4
    n = int(time/dt)#10**5
    planet_masses = system.masses
    star_mass = system.star_mass
    r = np.zeros((n,2))
    v = np.zeros_like(r)
    a = np.zeros_like(r)
    r[0] = position
    v[0] = velocity
    G = 4*np.pi
    dist_p6 = np.zeros(n)
    #start is timestep when rocket is launched
    start = 35024
    for i in trange(n-1):#calculate numerically rocket position with gravity from all planets
        for j in range(planet_positions[0].shape[0]):    
            R = r[i] - planet_positions[start+i,j]
            R_norm = np.linalg.norm(R)
            a[i] -= G*planet_masses[j]*R/R_norm**3


        a[i] -= G*star_mass*r[i]/np.linalg.norm(r[i])**3

        v[i+1] = v[i] + a[i]*dt

        r[i+1] = r[i] + v[i+1]*dt
        #finds distance between rocket and planet6
        dist_p6[i+1] = np.linalg.norm(planet_positions[start+i,6]-r[i+1])
   
    #find index of when distnce was smallest
    dist_p6[0] = dist_p6[1]
    min_dist = np.min(dist_p6)
    idx = np.where(min_dist == dist_p6)[0]


    plt.scatter(planet_positions[start+idx,6,0],planet_positions[start+idx,6,1],color='green')
    print(planet_positions[start+idx,6])

    #l is when we are close enough to begin landing sequence
    plt.scatter(r[idx,0],r[idx,1],color='blue')
    l = np.linalg.norm(r[idx])*np.sqrt(planet_masses[6]/(10*star_mass))
    print('l=',l)
    print('dist=',dist_p6[idx])
    print('idx=',idx)
    plt.plot(r[:,0],r[:,1],label='rocket_trajectory')





    end = start + int(1/dt*time)#plot planet trajectory of planet 0 and planet 6 for as long as the rocket is moving between them
    plt.plot(planet_positions[start:end,0,0],planet_positions[start:end,0,1],label='p0')
    plt.plot(planet_positions[start:end,6,0],planet_positions[start:end,6,1],label='p6')
    plt.axis('equal')
    plt.legend()
    plt.show()


planet_positions_at_t1 = planet_positions[35024]

r0 = np.array([-5.30521215 , 3.72464034])

v0 = np.array([-3.86666154 ,-2.14646642])

print(v0)


#find angle between planets for when to launch with hohmann transer angle
def hohmann_transfer_angle_to_timestep():
    r1 = system.semi_major_axes[0]
    r2 = system.semi_major_axes[6]

    angle = np.pi*(1-1/(2*np.sqrt(2)*(r1/r2+1)**(3/2)))
    print(angle)
    p0 = planet_positions[:,0]
    p6 = planet_positions[:,6]
    tol = 1e-3
    for i in range(86143):
        a = np.arctan2(p6[i,1]-p0[i,1],p6[i,0]-p0[i,0])
        if i == 15993:
            print(a*180/np.pi)
        if abs(a-angle) < tol:
            print('i=',i)
            return i
print(hohmann_transfer_angle_to_timestep())
print(utils.m_to_AU(5.7358e+09))
