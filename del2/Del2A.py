#EGEN KODE
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from numba import njit
from ast2000tools import utils
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)

Aphelion = system.aphelion_angles
Omega = system.initial_orbital_angles
A = system.semi_major_axes
E = system.eccentricities
n = 10000
theta = np.linspace(0,2*np.pi,n)


#LIST OF INITAL ORBIT ANGLES
''' throughout the project i use einsum alot. I mostly just use it to switch the indexing of the arrays and matrixes.
Thats why i usually never comment on it'''
f = np.einsum('ij->ji',np.linspace(system.initial_orbital_angles, 2*np.pi + system.initial_orbital_angles))

farger = ['brown','blue','pink','red','yellow','green','black']

#ANALYTIC SOLUTION
def R(a,e,f):
    return a*(1-e**2)/(1-e*np.cos(f))

#CALCULATE ANALYTIC ORBIT FOR ALL PLANETS
for planet_idx in range(system.number_of_planets):
    omega = Omega[planet_idx]
    a = A[planet_idx]
    e = E[planet_idx]
    aph = Aphelion[planet_idx]
    x = R(a,e,f[planet_idx])*np.cos(f[planet_idx]+aph)
    y = R(a,e,f[planet_idx])*np.sin(f[planet_idx]+aph)

 
    plt.plot(x,y,label=f'{planet_idx} a ',color=f'{farger[planet_idx]}')


# plt.scatter(0,0)
# plt.legend()
# plt.xlabel('$r_x$ in AU')
# plt.ylabel('$r_y$ in AU')
# plt.title('Analytical')
# plt.show()


n = 2*10**6
v = np.zeros((n,7,2))
r = np.zeros((n,7,2))
dt = 1e-4

v[0] = np.einsum('ij->ji',system.initial_velocities)
r[0] = np.einsum('ij->ji',system.initial_positions)

yr = 0
star_mass = system.star_mass
#INTEGRATION LOOP
@njit
def leapfrog(v,r,i,dt,yr,star_mass):
    G = 4*np.pi**2
    for j in range(7):
        r_norm = np.linalg.norm(r[i,j])

        a = -G*star_mass*(r[i,j]/r_norm**3)
        vi05 = v[i,j] + a*dt/2
        r[i+1,j] = r[i,j] + vi05*dt

        r_norm = np.linalg.norm(r[i+1,j])

        a =  -G*star_mass*(r[i+1,j]/r_norm**3)
        v[i+1,j] = vi05 + a*dt/2

        if j == 0:#COUNT YEARS OF PLANET 0
            if r[i+1,j,1] > 0 and r[i,j,1] < 0:
                yr += 1

    return v,r,yr


for i in trange(n-1):
    v,r,yr = leapfrog(v,r,i,dt,yr,star_mass)
    if yr >= 22:#BREAK AFTER 22YEARS
        break

print(yr)
print(i)        
r[i:] = r[i]#chop out useless zeros

# for f in range(len(farger)):
#     plt.plot(r[:,f,0],r[:,f,1],label=f'{f}',color=farger[f])

# plt.legend()
# plt.scatter(0,0)
# plt.xlabel('$r_x$ in AU')
# plt.ylabel('$r_y$ in AU')
# plt.title('Numeric')
# plt.show()
np.save('velocities.npy',v)

mission.verify_planet_positions(dt*n, np.einsum('ijk->kji',r))
#system.generate_orbit_video(np.linspace(0,dt*n,n),np.einsum('ijk->kji',r))
# 21 år
# 273187 tidssteg: 273187/21 > 10000

