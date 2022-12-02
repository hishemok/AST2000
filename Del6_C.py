#EGEN KODE
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from numba import njit
from scipy.interpolate import interp1d
from ast2000tools import utils
from ast2000tools import constants
#utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)

G = constants.G
planet_mass = system.masses[6]*constants.m_sun
m_p = constants.m_p
rho_surf = system.atmospheric_densities[6] 
planet_radius = system.radii[6]*1e3
gamma = 1.4
mu = 29
k_B = constants.k_B


@njit
def loop(r,T,rho,i,mu,gamma,m_p,pmass,prad,G):
    #integrate temperature and density by integrating the analytic solution
    grav = G*pmass/(prad+r[i])**2
    if T[i] > T[0]/2:
        dT_dr = -(gamma-1)/gamma*mu*m_p*grav/k_B
        T[i+1] = T[i] + dT_dr
        rho[i+1] = rho[i]  -rho[i]*dT_dr/T[i] -rho[i]*grav*mu*m_p/(T[i]*k_B)
    else:
        T[i+1] = T[0]/2
        rho[i+1] = rho[i] - rho[i]*grav*mu*m_p/(T[i+1]*k_B)
    return T,rho

def atmosphere():
    max_height = 300*1e3#top of atmosphere
    #dr = 1
    n = int(max_height)
    r = np.linspace(0,max_height,n+1)
    T = np.zeros(n+1)
    T[0] = 334.708
    rho = np.zeros(n+1)
    rho[0] = rho_surf
    for i in trange(n):
        T,rho = loop(r,T,rho,i,mu,gamma,m_p,planet_mass,planet_radius,G)
    return T,rho,r

temperature,density,r = atmosphere()
f1 = interp1d(r,temperature,kind='quadratic')
f2 = interp1d(r,density,kind='quadratic')
np.save('densities.npy',f2(r))
plt.plot(r/10**3,f1(r)/temperature[0])
plt.xlabel('r [km]')
plt.ylabel('Temp T/T0 [K]')
plt.show()
plt.plot(r/10**3,f2(r)/density[0])
plt.xlabel('r [km]')
plt.ylabel('Density rho/rho0 [kg/m^3]')
plt.show()
