import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from numba import njit

n = 10000
theta = np.linspace(0,2*np.pi,n)


from ast2000tools import utils
utils.check_for_newer_version()

seed = utils.get_seed('<Claudieg>')#Claudieg

from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)
Aphelion = system.aphelion_angles
Omega = system.initial_orbital_angles
A = system.semi_major_axes
E = system.eccentricities

f = np.einsum('ij->ji',np.linspace(system.initial_orbital_angles, 2*np.pi + system.initial_orbital_angles))


farger = ['brown','blue','pink','red','yellow','green','black']


def R(a,e,aph,f):
    return a*(1-e**2)/(1-e*np.cos(np.pi-aph+f))#

for planet_idx in range(system.number_of_planets):
    omega = Omega[planet_idx]
    a = A[planet_idx]
    e = E[planet_idx]
    aph = Aphelion[planet_idx]
    x = R(a,e,aph,f[planet_idx])*np.cos(f[planet_idx])
    y = R(a,e,aph,f[planet_idx])*np.sin(f[planet_idx])

 
    plt.plot(x,y,label=f'{planet_idx} a',color=f'{farger[planet_idx]}')


# plt.legend()
# plt.scatter(0,0)
# plt.show()


n = 10**6
v = np.zeros((n,7,2))
r = np.zeros((n,7,2))
dt = 6.5*1e-6

v[0] = np.einsum('ij->ji',system.initial_velocities)
r[0] = np.einsum('ij->ji',system.initial_positions)

yr = 0
star_mass = system.star_mass
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

        if j == 0:
            if r[i+1,j,1] > 0 and r[i,j,1] < 0:
                yr += 1

    return v,r,yr


for i in trange(n-1):
    v,r,yr = leapfrog(v,r,i,dt,yr,star_mass)
    if yr >= 22:
        break

print(yr)
print(i)        
r[i:] = r[i]
for f in range(len(farger)):
    plt.plot(r[:,f,0],r[:,f,1],label=f'{f}',color=farger[f])

plt.legend()
plt.scatter(0,0)
plt.show()


system.verify_planet_positions(dt*n,np.einsum('ijk->kji',r))
system.generate_orbit_video(np.linspace(0,dt*n,n),np.einsum('ijk->kji',r))
# 22 år
# 273187 tidssteg: 273187/21 > 10000

