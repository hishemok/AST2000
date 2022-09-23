from lib2to3.pgen2.token import RPAR
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from numba import njit


from ast2000tools import utils
utils.check_for_newer_version()
seed = utils.get_seed('<Claudieg>')
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)



n = int(1.5*10**4)
r = np.zeros((n,2,2))
v = np.zeros(r.shape)
dt = 1e-4

r[0,0] = np.array([0,0])
r[0,1] = np.einsum('ij->ji',system.initial_positions)[0]
v[0,1] = np.einsum('ij->ji',system.initial_velocities)[0]


print(np.linalg.norm(np.einsum('ij->ji',system.initial_positions),axis=1))


m_sun = system.star_mass
m_planet = system.masses[0]
G = 4*np.pi**2




#P_p = -P_s
v[0,0] = -v[0,1]*m_planet/m_sun

for i in trange(n-1):
    r_norm = np.linalg.norm(r[i,1]-r[i,0])
    F = -G*m_sun*m_planet*(r[i,1]-r[i,0])/r_norm**3

    vi05_s = v[i,0] - F*dt/2/m_sun
    vi05_p = v[i,1] + F*dt/2/m_planet

    r[i+1,0] = r[i,0] + vi05_s*dt
    r[i+1,1] = r[i,1] + vi05_p*dt

    r_norm = np.linalg.norm(r[i+1,1]-r[i+1,0])
    F = -G*m_sun*m_planet*(r[i+1,1]-r[i+1,0])/r_norm**3

    v[i+1,0] = vi05_s - F*dt/2/m_sun
    v[i+1,1] = vi05_p + F*dt/2/m_planet



radius_planet = utils.m_to_AU(system.radii[0]*10**3)
radius_star = utils.m_to_AU(system.star_radius*10**3)

#To simplify assume they are squared
A_planet = np.pi*radius_planet**2
A_star = np.pi*radius_star**2
obs_dist = 10#AU
d1 = radius_planet
r1 = obs_dist-max(np.einsum('ij->i',r[:,1]))
r2 = obs_dist
scale = r2/r1
#print(scale)




def shadowline():
    angle_s = radius_star/obs_dist
    dp = angle_s*r1
    n = 1000
    r = np.linspace(-2*dp,2*dp,n)
    diameter_planet = radius_planet*2
    #planet visible within  -dp < r < dp 
    shadow = np.zeros(n)
    for i in range(len(r)):
        if r[i]-radius_planet > -dp and r[i]+radius_planet < dp:
            shadow[i] = 1
        elif r[i]+radius_planet < -dp or r[i]-radius_planet > dp:
            shadow[i] = 0
        elif r[i] + radius_planet > -dp and r[i] - radius_planet < -dp:
            L = r[i]+radius_planet+dp
            shadow[i] = 1-(diameter_planet-L)/diameter_planet
        elif r[i] + radius_planet > dp and r[i] - radius_planet < dp:
            L = r[i]+radius_planet-dp
            shadow[i] = (diameter_planet-L)/diameter_planet
            
    return shadow




shadow = shadowline()

def f(line):
    return (A_star-scale*A_planet*line)/A_star

n_max = np.amax(f(shadow))
standard_deviation = 1e-4
noice =np.random.normal(-standard_deviation,standard_deviation,1000)

flux = f(shadow)+noice

plt.plot(flux,label='with_noice')
# plt.plot(f(shadow),label='no_noice')
plt.legend()
plt.show()

