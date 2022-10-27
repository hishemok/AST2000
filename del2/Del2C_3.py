import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from numba import njit


from ast2000tools import utils
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)
from ast2000tools import constants



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
# print(max(F))
# print(F)
#largest gravitational force from last planet

n = int(2*10**4)
r = np.zeros((n,2,2))
v = np.zeros(r.shape)
dt = 1e-4

r[0,0] = np.array([0,0])
r[0,1] = np.einsum('ij->ji',system.initial_positions)[-1]
v[0,1] = np.einsum('ij->ji',system.initial_velocities)[-1]

m_sun = system.star_mass
m_planet = system.masses[-1]
G = 4*np.pi**2

#conserving momentum
v[0,0] = -v[0,1]*m_planet/m_sun

Ep = np.zeros(n-1)
Ek = np.zeros_like(Ep)



def E_tot(v1,v2,r1,r2):
    M = system.star_mass
    m = system.masses[-1]
    m_red= M*m/(M+m)
    v = np.linalg.norm(v1)+np.linalg.norm(v2)
    r = np.linalg.norm(r1-r2)

    G = 4*np.pi**2
    return 1/2*m_red*v**2, - G*M*m/r



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


    


plt.plot(v[:,0,0],v[:,0,1])
plt.plot(v[:,1,0],v[:,1,1])
plt.axis('equal')
plt.show()

angle = np.pi/3
v_pec = 1e-5
v_rad = np.einsum('ij->i',(v[:,0]+v_pec)*np.sin(angle))
v_max = np.amax(v_rad)
noise = np.random.normal(0,0.2*v_max,n)#Draw random samples from a normal (Gaussian) distribution.

plt.plot(np.linspace(0,n*dt,n),v_rad+noise,label='v_rad')
plt.xlabel('t in yrs')
plt.ylabel('$v_rad$')
plt.title('Radial velocity')
plt.legend()
plt.show()


