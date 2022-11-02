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




'''
def rocket_trajectory(planet_positions_at_t,position,velocity,time):
    dt = 1e-4
    n = int(time/dt)#10**5
    planet_masses = system.masses
    star_mass = system.star_mass
    r_planets = planet_positions_at_t

    # dt = time/n

    r = np.zeros((n,2))
    v = np.zeros_like(r)
    a = np.zeros_like(r)
    r[0] = position
    v[0] = velocity
    G = 4*np.pi
    dist_p6 = np.zeros(n)
    #start is timestep when rocket is launched
    start = 15993
    print(planet_positions[15993+69280,6])

    for i in trange(n-1):
        for j in range(r_planets.shape[0]):
            
            R = r[i] - r_planets[j]
            R_norm = np.linalg.norm(R)
            a[i] -= G*planet_masses[j]*R/R_norm**3


        a[i] -= G*star_mass*r[i]/np.linalg.norm(r[i])**3
        v[i+1] = v[i] + a[i]*dt
        #without any speedboost this timestep was the closest between whem
        if i == 68607:
            v[i+1] += np.array([-.6,-1.05])
        #second correction
        if i == 68699:
            v[i+1] += np.array([.4,-1])
        r[i+1] = r[i] + v[i+1]*dt
        #calculates distance between rocket and planet6
        dist_p6[i+1] = np.linalg.norm(planet_positions[start+i,6]-r[i+1])
   
    #find index of when distnce was smallest
    dist_p6[0] = dist_p6[1]
    min_dist = np.min(dist_p6)
    idx = np.where(min_dist == dist_p6)[0]
    #at idx = 69280, distance was d = 0.00016068 < l=0.001857530660315445
    #where l is required to land at given planet
    #remember idx+start is the true timestep

    plt.scatter(planet_positions[start+idx,6,0],planet_positions[start+idx,6,1])
    #[-12.7432317] [5.38150351] planeten
    plt.scatter(r[idx,0],r[idx,1])
    l = np.linalg.norm(r[idx])*np.sqrt(planet_masses[6]/(10*star_mass))
    print(l)
    print(dist_p6[idx])
    print(idx)



    plt.plot(r[:,0],r[:,1],label='rocket_trajectory')
    
    end = start + int(1/dt*time)
    plt.plot(planet_positions[start:end,0,0],planet_positions[start:end,0,1],label='p0')
    plt.plot(planet_positions[start:end,6,0],planet_positions[start:end,6,1],label='p6')
    plt.axis('equal')
    plt.legend()
    plt.show()
'''
    

def rocket_trajectory_2(position,velocity,time):
    dt = 1e-4
    n = int(time/dt)#10**5
    planet_masses = system.masses
    star_mass = system.star_mass

    # dt = time/n

    r = np.zeros((n,2))
    v = np.zeros_like(r)
    a = np.zeros_like(r)
    r[0] = position
    v[0] = velocity
    G = 4*np.pi
    dist_p6 = np.zeros(n)
    #start is timestep when rocket is launched
    start = 35024
    for i in trange(n-1):
        for j in range(planet_positions[0].shape[0]):    
            R = r[i] - planet_positions[start+i,j]
            R_norm = np.linalg.norm(R)
            a[i] -= G*planet_masses[j]*R/R_norm**3


        a[i] -= G*star_mass*r[i]/np.linalg.norm(r[i])**3

        
        # if i < 15:
        #     print(a[i])

        v[i+1] = v[i] + a[i]*dt

        r[i+1] = r[i] + v[i+1]*dt

        dist_p6[i+1] = np.linalg.norm(planet_positions[start+i,6]-r[i+1])
   
    #find index of when distnce was smallest
    dist_p6[0] = dist_p6[1]
    min_dist = np.min(dist_p6)
    idx = np.where(min_dist == dist_p6)[0]


    plt.scatter(planet_positions[start+idx,6,0],planet_positions[start+idx,6,1],color='green')
    print(planet_positions[start+idx,6])


    plt.scatter(r[idx,0],r[idx,1],color='blue')
    l = np.linalg.norm(r[idx])*np.sqrt(planet_masses[6]/(10*star_mass))
    print('l=',l)
    print('dist=',dist_p6[idx])
    print('idx=',idx)
    plt.plot(r[:,0],r[:,1],label='rocket_trajectory')





    end = start + int(1/dt*time)
    plt.plot(planet_positions[start:end,0,0],planet_positions[start:end,0,1],label='p0')
    plt.plot(planet_positions[start:end,6,0],planet_positions[start:end,6,1],label='p6')
    plt.axis('equal')
    plt.legend()
    plt.show()

# r0 = np.array([2.86591894 ,6.01543821])
# v0 = np.array([-2.1817216 ,2.32774964])*0.976


planet_positions_at_t1 = planet_positions[35024]

r0 = np.array([-5.30521215 , 3.72464034])

v0 = np.array([-3.86666154 ,-2.14646642])

print(v0)


# rocket_trajectory_2(r0,v0,10)


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
print([0.165-1.74254273e-02-3.96789707e-03,-0.05+1.88926670e-05+7.11948886e-05])