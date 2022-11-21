#EGEN KODE
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools import utils
from ast2000tools import constants
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)
from ast2000tools import star_population

sigma = constants.sigma
c = constants.c
k_B = constants.k_B
G = constants.G

#def sun info
R_sun = constants.R_sun
L_sun = constants.L_sun
M_sun = constants.m_sun

#our star
star_radius = system.star_radius*1e3
star_temperature = system.star_temperature
star_mass = system.star_mass

#A1
def HR_diagram():
    stars = star_population.StarPopulation()
    T = stars.temperatures # [K]
    L = stars.luminosities # [L_sun]
    r = stars.radii        # [R_sun]

    c = stars.colors
    s = np.maximum(1e3*(r - r.min())/(r.max() - r.min()), 1.0) # Make point areas proportional to star radii

    fig, ax = plt.subplots()
    
    Luminocity = (4*np.pi*star_radius**2*sigma*star_temperature**4)/L_sun

    ax.scatter(star_temperature,Luminocity,color='red',label='Our star')
    ax.scatter(T, L, c=c, s=s, alpha=0.8, edgecolor='k', linewidth=0.05)

    ax.set_xlabel('Temperature [K]')
    ax.invert_xaxis()
    ax.set_xscale('log')
    ax.set_xticks([35000, 18000, 10000, 6000, 4000, 3000])
    ax.set_xticklabels(list(map(str, ax.get_xticks())))
    ax.set_xlim(40000, 2000)
    ax.minorticks_off()

    ax.set_ylabel(r'Luminosity [$L_\odot$]')
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1e6)
    plt.savefig('HR_diagram.png')

    plt.legend()
    plt.show()
    return Luminocity

Luminocity = HR_diagram()

#A2
T_mainsequence = 0.1*star_mass*M_sun*c**2*0.007/(L_sun*Luminocity)
print(T_mainsequence)
'''
1.4498989360959194e+16
'''

#A3
print(Luminocity*L_sun/(star_mass*M_sun)**4)
print(L_sun/M_sun**4)
'''
9.452541315903319e-96
2.4484487919948118e-95
pretty good
'''

print(star_temperature*star_radius/(star_mass*M_sun))
print(5778*R_sun/M_sun)
'''
2.666349749776749e-18
2.0215259233247465e-18
well behaved :)
'''

#B1
H_m = 1.6735575*1e-27
HE_m = 6.6464731*1e-27
T = 10 #K
mu = (H_m*0.75+HE_m*0.25)/H_m
#GMC Radius
Rad = G*star_mass*M_sun*mu*H_m/(5*k_B*T)
print(Rad)
print(utils.m_to_AU(Rad))

def HR_diagram2():
    print('HR2')
    H_m = 1.6735575*1e-27
    HE_m = 6.6464731*1e-27
    T = 10 #K
    mu = (H_m*0.75+HE_m*0.25)/H_m
    #GMC Radius
    Rad = G*star_mass*M_sun*mu*H_m/(5*k_B*T)
    #GMC Luminocity
    Luminocity_GMC = (4*np.pi*Rad**2*sigma*T**4)/L_sun
    print(Luminocity_GMC)


    stars = star_population.StarPopulation()
    T = stars.temperatures # [K]
    L = stars.luminosities # [L_sun]
    r = stars.radii        # [R_sun]

    c = stars.colors
    s = np.maximum(1e3*(r - r.min())/(r.max() - r.min()), 1.0) # Make point areas proportional to star radii

    fig, ax = plt.subplots()
    

    ax.scatter(T, L, c=c, s=s, alpha=0.8, edgecolor='k', linewidth=0.05)

    Luminocity = (4*np.pi*star_radius**2*sigma*star_temperature**4)/L_sun
    ax.scatter(star_temperature,Luminocity,color='red',label='Our star')
    print(Luminocity)

    ax.scatter(10,Luminocity_GMC,color='g',label='GMC')

    ax.set_xlabel('Temperature [K]')
    ax.invert_xaxis()
    ax.set_xscale('log')
    ax.set_xticks([35000, 10000, 1000, 500 , 9])
    ax.set_xticklabels(list(map(str, ax.get_xticks())))
    ax.set_xlim(40000, 9)
    ax.minorticks_off()

    ax.set_ylabel(r'Luminosity [$L_\odot$]')
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1e6)
    plt.savefig('HR_diagram.png')

    plt.legend()
    plt.show()

HR_diagram2()


