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


r = np.einsum('ij->ji',system.initial_positions)
m = system.masses

#FUNCTION CALCULATES GFORCE OF GIVEN PLANET
def Gforce(r_vec,m_planet):
    r_norm = np.linalg.norm(r_vec)
    m_sun = 0.316341
    G = 4*np.pi**2
    return -G*m_planet*m_sun/r_norm**2

F = np.zeros(7)#F IS LIST OF GFORCE
for planet in range(system.number_of_planets):
    F[planet] = Gforce(r[planet],m[planet])

F = np.abs(F)#MAKE THEM POSITIVE AND FIND GREATEST
# print(max(F))
# print(F)
#størst Fg fra siste planet

n = 2*10**4
r = np.zeros((n,2,2))
v = np.zeros(r.shape)
dt = 1e-4

r[0,0] = np.array([0,0])
r[0,1] = np.einsum('ij->ji',system.initial_positions)[-1]
v[0,1] = np.einsum('ij->ji',system.initial_velocities)[-1]

m_sun = system.star_mass
m_planet = system.masses[-1]
G = 4*np.pi**2

#P_p = -P_s
v[0,0] = -v[0,1]*m_planet/m_sun

E = np.zeros(n-1)
Ep = np.zeros_like(E)
Ek = np.zeros_like(E)


def E_tot(v1,v2,r1,r2):#CALCULATES TOTAL ENERGY OF PLANET IN ORBIT
    M = system.star_mass
    m = system.masses[-1]
    m_red= M*m/(M+m)
    v = np.linalg.norm(v1+v2)
    r = np.linalg.norm(r1+r2)
    G = 4*np.pi**2
    return 1/2*m_red*v**2, - G*M*m/r




for i in trange(n-1):# INTEGRATION LOOP 
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

    
    Ek[i],Ep[i] = E_tot(v[i+1,1],v[i+1,0],r[i+1,1],r[i+1,0])
    


plt.plot(r[:,0,0],r[:,0,1])
plt.plot(r[:,1,0],r[:,1,1])
plt.xlabel('x/AU')
plt.ylabel('y/AU')
plt.show()
plt.plot(Ek,label='Kinetic')
plt.plot(Ep,label='Potential')
plt.plot(Ek+Ep,label='Total')
plt.xlabel('time')
plt.ylabel('Energy')
#PLOT POTENTIAL AND KINETIC AND TOTAL ENERGY

plt.legend()
plt.show()

