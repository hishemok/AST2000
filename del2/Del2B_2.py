import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from numba import njit


from ast2000tools import utils
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)

n = 100
theta = np.linspace(0,2*np.pi,n)

def R(a,e,aph,f):
    return a*(1-e**2)/(1-e*np.cos(np.pi-aph+f))

def Kepler_newton(index):
    a = system.semi_major_axes[index]
    e = system.eccentricities[index]
    m = system.masses[index]
    star_mass = 0.316341
    b = a*np.sqrt(1-e**2)
    h = np.sqrt(a*(1-e**2)*m)
    P = 2*np.pi*a*b/h

    G = 4*np.pi**2

    Kepler = np.array([P**2,a**3])
    Newton = np.array([P**2,4*np.pi**2*a**3/(G*(m+star_mass))])
    return Kepler,Newton


for i in range(7):
    Kepler,Newton = Kepler_newton(i)
    #print(f'P^2 {Kepler[0]},     Keplers a^3 {Kepler[1]},    Newtons {Newton[1]}')













n = 2*10**6
v = np.zeros((n,7,2))
r = np.zeros((n,7,2))
dt = 1e-4

v[0] = np.einsum('ij->ji',system.initial_velocities)
r[0] = np.einsum('ij->ji',system.initial_positions)

A = system.semi_major_axes

P = np.zeros(7)
P_ = np.zeros(7)

yr = 0
star_mass = system.star_mass
@njit
def leapfrog(v,r,i,dt,yr,P,P_,star_mass):
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


        if P[j] == 0 and P_[j] == 0:
            if r[i+1,j,1] > 0 and r[i,j,1] < 0:
                P[j] = dt*i
                P_[j] = 1
        elif P[j] != 0 and P_[j] == 1:
            if r[i+1,j,1] > 0 and r[i,j,1] < 0:
                P[j] -= dt*i
                P_[j] = 2  



    return v,r,yr,P,P_


for i in trange(n-1):
    v,r,yr,P,P_ = leapfrog(v,r,i,dt,yr,P,P_,star_mass)
    if yr >= 42:
        break



for planet in range(system.number_of_planets):
    p = P[planet]**2
    k = A[planet]**3
    n = A[planet]**3/(system.masses[planet]+star_mass)
    print(f'Planet {planet}: P^2= {p:.7f}   K:{k:.7f}    N:{n:.7f}')



# r[i:] = r[i]
# for planet_idx in range(system.number_of_planets):
#     plt.plot(r[:,planet_idx,0],r[:,planet_idx,1],label=f'{planet_idx}')

# # plt.legend()
# plt.scatter(0,0)
# plt.show()


# 21 år
# 273187 tidssteg: 273187/21 > 10000

