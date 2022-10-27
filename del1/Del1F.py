from ast2000tools import utils
utils.check_for_newer_version()

from ast2000tools import constants as const
print('One astronomical unit is defined as', const.AU, 'meters.')

seed = utils.get_seed('Claudieg')

from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)

print('My system has a {:g} solar mass star with a radius of {:g} kilometers.'
      .format(system.star_mass, system.star_radius))

for planet_idx in range(system.number_of_planets):
    print('Planet {:d} is a {} planet with a semi-major axis of {:g} AU.'
          .format(planet_idx, system.types[planet_idx], system.semi_major_axes[planet_idx]))


for planet_idx in range(system.number_of_planets):
    print('Planet {:d}  has orbital angle {:g} rad. Eccentricities'
          .format(planet_idx,  system.initial_orbital_angles[planet_idx]), system.eccentricities[planet_idx])


from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)

# home_planet_idx = 0 # The home planet always has index 0
# print('My mission starts on planet {:d}, which has a radius of {:g} kilometers.'
#       .format(home_planet_idx, mission.system.radii[home_planet_idx]))

# print('My spacecraft has a mass of {:g} kg and a cross-sectional area of {:g} m^2.'
#       .format(mission.spacecraft_mass, mission.spacecraft_area))

# i = 0
# for planet_idx in range(system.number_of_planets):
#     print(f'Planet with index {i} has a mass of: {system.masses[planet_idx]*1.989e30}')
#     i+=1

system.print_info()

#print(utils.AU_pr_yr_to_m_pr_s(11.51))
print(utils.AU_to_m(8.20895e-06))


