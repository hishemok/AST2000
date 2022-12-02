#EGEN KODE
#%%
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

'''B'''
#1
G = constants.G
sigma = constants.sigma
m_p = constants.m_p

star_mass = system.star_mass
M_sun = constants.m_sun
M_ch = M_sun*1.4
R_sun = constants.R_sun


M_wd = star_mass*M_ch/(M_sun*8) * M_sun#last m_sun to turn from solar mass to kg



Z,A = 1,2 #Z/A = 0.5
h =  6.62607015*10**(-34)
m_e = 9.109*10**(-31)

R_wd = (3/2*np.pi)**(4/3) * h**2/(20*m_e*G) * (Z/(A*m_p))**(5/3) * M_wd**(-1/3)/R_sun# in sun radii
print('radius white dwarf:  ',R_wd)
#0.0495  sun radii

#2
R_wd *= R_sun #back to meters
V = 4/3*R_wd**3
density = M_wd/R_wd #kg/m^3
density_litre = density/1000 #litre from m^3s
print('White dwarf density      ',density_litre)
#3.9133532743332135e+19 damn
#10^10 ganger for stor så vi må gjøre noe med den saken

surface_grav_acc = np.sqrt(G*M_wd/R_wd)
print('surface grav acc     ',surface_grav_acc)

'''----------------------'''

star_mass = star_mass
star_temperature = system.star_temperature
star_radius = system.star_radius*1e3
L_star = sigma*star_temperature**4*4*np.pi*star_temperature**2

M_sun = constants.m_sun
L_sun = constants.L_sun


stars = star_population.StarPopulation()
T = stars.temperatures # [K]
L = stars.luminosities # [L_sun]
r = stars.radii        # [R_sun]

c = stars.colors
s = np.maximum(1e3*(r - r.min())/(r.max() - r.min()), 1.0) # Make point areas proportional to star radii
s_sun = np.maximum(1e3*(star_radius/R_sun-r.min())/(r.max()-r.min()),1.0)

fig, ax = plt.subplots()

Luminocity = (4*np.pi*star_radius**2*sigma*star_temperature**4)/L_sun

ax.scatter(T, L, c=c, s=s, alpha=0.8, edgecolor='k', linewidth=0.05)
ax.scatter(star_temperature,Luminocity,color='red',label='Our star')

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


'''------------------------------'''
''' Plot Sub Giant in HR diagram '''
fig,ax = plt.subplots()
ax.scatter(T, L, c=c, s=s, alpha=0.8, edgecolor='k', linewidth=0.05)
ax.scatter(star_temperature,Luminocity,color='red',label='Our star')
ax.scatter(4000, L_star/L_sun+140,alpha=0.8, s=s_sun*50, color = 'blue',label='Sub giant')
ax.legend()
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
plt.show()

''' Plot Red Giant in HR diagram '''
fig,ax = plt.subplots()
ax.scatter(T, L, c=c, s=s, alpha=0.8, edgecolor='k', linewidth=0.05)
ax.scatter(star_temperature,Luminocity,color='red',label='Our star')
ax.scatter(2600, 1e3,s=s_sun*100,alpha=0.8,color = 'blue',label='red giant')
ax.legend()
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
plt.show()

''' Plot Super Giant in HR diagram '''
fig,ax = plt.subplots()
ax.scatter(T, L, c=c, s=s, alpha=0.8, edgecolor='k', linewidth=0.05)
ax.scatter(star_temperature,Luminocity,color='red',label='Our star')
ax.scatter(2600, 1e5,s=s_sun*300,alpha=0.8,color = 'blue',label='Super giant')
ax.legend()
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
plt.show()


''' Plot White Dwarf in HR diagram '''
fig,ax = plt.subplots()
ax.scatter(T, L, c=c, s=s, alpha=0.8, edgecolor='k', linewidth=0.05)
ax.scatter(star_temperature,Luminocity,color='red',label='Our star')
ax.scatter(30000, 10**(-1.4),s=s_sun*5,color = 'blue',label='white dwarf')
ax.legend()
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
plt.show()



