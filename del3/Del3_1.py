#EGEN KODE
import numpy as np
from ast2000tools import constants
from ast2000tools import utils
#utils.check_for_newer_version()

seed = utils.get_seed('Claudieg')#Claudieg

from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)

sigma = constants.sigma
T = system.star_temperature

r = utils.AU_to_m(np.linalg.norm(np.einsum('ij->ji',system.initial_positions)[0]))
rs = system.star_radius*10**3 #m
A = 4*np.pi*rs**2
dA = A/(4*np.pi*r**2)


F = (rs/r)**2*sigma*T**4
W = 40/0.12
Ap = 2*np.pi*(system.radii[0]*10**3)**2
#print(F*Ap)
print(W/F)

def flux_on_planet(r,R):
    F = (rs/r)**2*sigma*T**4
    A = 2*np.pi*R**2
    return F*A

#test planet 0
#print(flux_on_planet(r,system.radii[0]*10**3))
#same as above

def planet_temperature(flux):
    return (flux/(sigma*2*np.pi*R**2))**(1/4)


habitable = []
for planet_idx in range(system.number_of_planets):
    r = utils.AU_to_m(np.linalg.norm(np.einsum('ij->ji',system.initial_positions)[planet_idx]))
    R = system.radii[0]*10**3
    flux = flux_on_planet(r,R)
    print(f'Planet {planet_idx}: T {planet_temperature(flux):.3f} K')
    if planet_temperature(flux) > 245 and planet_temperature(flux) < 410:
        # print(f'Planet {planet_idx} is habitable ')
        habitable.append(planet_idx)
print('')
for i in habitable:
    print(f'Planet {i} is habitable')

    density = (system.masses[i]*constants.m_sun)/((4/3)*np.pi*(system.radii[i]*10**3)**3) #mass/volume
    print(f"density of planet {i} is: {density}kg/m^3")


print(f"Temperature of sun: {T}")
print(f"radius of sun: {rs}")
print(f"Surface area of sun: {A}")
print(f"luminosity of sun: {A*sigma*T**4}")
print(f"Surface flux of sun: {sigma*T**4}")

print(f"distance of planet 4 from sun is: {utils.AU_to_km(np.linalg.norm(np.einsum('ij->ji',system.initial_positions)[4]))}km")
print(f"distance of planet 6 from sun is: {utils.AU_to_km(np.linalg.norm(np.einsum('ij->ji',system.initial_positions)[6]))}km")
print(f"distance of planet 0 from sun is: {utils.AU_to_km(np.linalg.norm(np.einsum('ij->ji',system.initial_positions)[0]))}km")

print(f"distance of planet 2 from sun is: {utils.AU_to_km(np.linalg.norm(np.einsum('ij->ji',system.initial_positions)[2]))}km")
print(f"radius 0: {system.radii[0]*10**3} radius 4: {system.radii[4]*10**3} radius 6: {system.radii[6]*10**3} radius 2: {system.radii[2]*10**3}")