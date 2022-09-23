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


r = np.einsum('ij->ji',system.initial_positions)
m = system.masses

def Gforce(r_vec,m_planet):
    r_norm = np.linalg.norm(r_vec)
    m_sun = 0.316341
    G = 4*np.pi**2
    return -G*m_planet*m_sun/r_norm**2

F = np.zeros(7)
for planet in range(system.number_of_planets):
    F[planet] = Gforce(r[planet],m[planet])

F = np.abs(F)

# print(F)
#Three planets whith greatest Gforce
#index(0,5,6)


# print(F)
#størst Fg fra siste planet

n = int(1.5*10**4)
r = np.zeros((n,4,2))
v = np.zeros(r.shape)
dt = 1e-4

r[0,0] = np.array([0,0])
r[0,1] = np.einsum('ij->ji',system.initial_positions)[0]
v[0,1] = np.einsum('ij->ji',system.initial_velocities)[0]
r[0,2] = np.einsum('ij->ji',system.initial_positions)[-2]
v[0,2] = np.einsum('ij->ji',system.initial_velocities)[-2]
r[0,3] = np.einsum('ij->ji',system.initial_positions)[-1]
v[0,3] = np.einsum('ij->ji',system.initial_velocities)[-1]

m_sun = system.star_mass
m_planet0 = system.masses[0]
m_planet1 = system.masses[-2]
m_planet2 = system.masses[-1]
G = 4*np.pi**2

planet_masses = [m_planet0,m_planet1,m_planet2]
#P_p = -P_s
v[0,0] = -(v[0,1]*m_planet0 +v[0,2]*m_planet0 + v[0,3]*m_planet2)/m_sun


a = np.zeros_like(r)
P = 0
vel = np.zeros(n-1)
for i in trange(n-1):
    r_norm1 = np.linalg.norm(r[i,1]-r[i,0])
    r_norm2 = np.linalg.norm(r[i,2]-r[i,0])
    r_norm3 = np.linalg.norm(r[i,3]-r[i,0])

    F1 = -G*m_sun*planet_masses[0]*(r[i,1]-r[i,0])/r_norm1**3
    F2 = -G*m_sun*planet_masses[1]*(r[i,2]-r[i,0])/r_norm2**3
    F3 = -G*m_sun*planet_masses[2]*(r[i,3]-r[i,0])/r_norm3**3

    v05_s = v[i,0] - F1*dt/2/m_sun - F2*dt/2/m_sun - F3*dt/2/m_sun
    v05_1 = v[i,1] + F1*dt/2/planet_masses[0]
    v05_2 = v[i,2] + F2*dt/2/planet_masses[1]
    v05_3 = v[i,3] + F3*dt/2/planet_masses[2]

    r[i+1,0] += r[i,0] + v05_s*dt
    r[i+1,1] = r[i,1] + v05_1*dt
    r[i+1,2] = r[i,2] + v05_2*dt
    r[i+1,3] = r[i,3] + v05_3*dt


    r_norm1 = np.linalg.norm(r[i+1,1]-r[i+1,0])
    r_norm2 = np.linalg.norm(r[i+1,2]-r[i+1,0])
    r_norm3 = np.linalg.norm(r[i+1,3]-r[i+1,0])

    F1 = -G*m_sun*planet_masses[0]*(r[i+1,1]-r[i+1,0])/r_norm1**3
    F2 = -G*m_sun*planet_masses[1]*(r[i+1,2]-r[i+1,0])/r_norm2**3
    F3 = -G*m_sun*planet_masses[2]*(r[i+1,3]-r[i+1,0])/r_norm3**3

    v[i+1,0] += v05_s - F1*dt/2/m_sun - F2*dt/2/m_sun - F3*dt/2/m_sun
    v[i+1,1] = v05_1 + F1*dt/2/planet_masses[0]
    v[i+1,2] = v05_2 + F2*dt/2/planet_masses[1]
    v[i+1,3] = v05_3 + F3*dt/2/planet_masses[2]







plt.plot(r[:,0,0],r[:,0,1])
plt.plot(r[:,1,0],r[:,1,1])
plt.plot(r[:,2,0],r[:,2,1])
plt.plot(r[:,3,0],r[:,3,1])

plt.show()


angle = np.pi/3
v_pec = np.linalg.norm(v[0,0])
v_rad = np.einsum('ij->i',(v[:,0]+v_pec)*np.sin(angle))
v_max = np.amax(v[:,0])
noice = np.einsum('ij->i',v[:,0]) + np.random.normal(-0.2*v_max,0.2*v_max,n)#Draw random samples from a normal (Gaussian) distribution.

plt.plot(v_rad+noice,label='v_rad')
plt.legend()
plt.show()

