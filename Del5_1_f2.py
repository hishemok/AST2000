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

planet_positions_at_t0 = np.einsum('ij->ji',system.initial_positions)

planet_trajectories = np.load('planet_trajectories.npz')
planet_positions = np.einsum('ijk->kji',planet_trajectories['planet_positions'])
planet_velocities = np.load('velocities.npy')

start = 76965#start timestep fra del5_1 hohmann angle

def rocket_trajectory(rocket_r0,rocket_v0,time):
    dt = 1e-5
    n = int(time/dt)#10**5
    planet_masses = system.masses
    star_mass = system.star_mass

    r = np.zeros((n,2))
    v = np.zeros_like(r)

    r[0] = rocket_r0
    v[0] = rocket_v0

    G = 4*np.pi**2
    dist_p6 = np.zeros(n)

    #im using njit because the rocket needs a very small dt for accurate measurements
    #there is a lot of mix between index i/10 and i, because i had to turn dt form 1e-4 to 1e-5
    @njit
    def loop(G,number_of_planets,r,v,i,dt,planet_positions):
        R = r[i]
        R_norm = np.linalg.norm(R)
        a = -G*star_mass*R/R_norm**3
        for j in range(number_of_planets):
            R = r[i] - planet_positions[start+int(i/10),j]
            R_norm = np.linalg.norm(R)
            a -= G*planet_masses[j]*R/R_norm**3
        
        v[i+1] = v[i] + a*dt
        if i == 1999900:
            v[i+1] += np.array([0,0.12])

        r[i+1] = r[i] + v[i+1]*dt
        distance = np.linalg.norm(planet_positions[start+int(i/10),6]-r[i+1])
        return r,v,distance

    for i in trange(n-1):
        r,v,dist_p6[i+1] = loop(G,system.number_of_planets,r,v,i,dt,planet_positions)
         

    dist_p6[0] = dist_p6[1]
    min_dist = np.min(dist_p6)
    idx = np.where(min_dist == dist_p6)[0]
    if idx.shape[0] > 1:
        idx = idx[1]
    
    plt.scatter(r[idx,0],r[idx,1],color='green')
    l = np.linalg.norm(r[idx])*np.sqrt(planet_masses[6]/(10*star_mass))
    print('l=',l)
    print('dist=',dist_p6[idx])
    print('idx=',idx)

    #with given speedboost at timestep = 1999900 dist<l at timestep 2077821
    plt.scatter(planet_positions[start+int(idx/10),6,0],planet_positions[start+int(idx/10),6,1],color='orange')
    
    plt.plot(planet_positions[start:start+int(i/10),0,0],planet_positions[start:start+int(i/10),0,1],label='p0')
    plt.plot(planet_positions[start:start+int(i/10),6,0],planet_positions[start:start+int(i/10),6,1],label='p6')
    plt.plot(r[:,0],r[:,1],label='rocket')

    plt.legend()
    plt.axis('equal')
    plt.show()




r0 = np.array([ 5.37926532 ,-4.02603141])
v0 = np.array([4.56978669 ,2.3374726 ])*1.163 + np.array([0.02,0])



rocket_trajectory(r0,v0,22)