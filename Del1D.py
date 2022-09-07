from ast2000tools import utils
utils.check_for_newer_version()

from ast2000tools import constants as const
print('One astronomical unit is defined as', const.AU, 'meters.')

seed = utils.get_seed('<Hishemok>')

from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)

print('My system has a {:g} solar mass star with a radius of {:g} kilometers.'
      .format(system.star_mass, system.star_radius))

for planet_idx in range(system.number_of_planets):
    print('Planet {:d} is a {} planet with a semi-major axis of {:g} AU.'
          .format(planet_idx, system.types[planet_idx], system.semi_major_axes[planet_idx]))

from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)

home_planet_idx = 0 # The home planet always has index 0
print('My mission starts on planet {:d}, which has a radius of {:g} kilometers.'
      .format(home_planet_idx, mission.system.radii[home_planet_idx]))

print('My spacecraft has a mass of {:g} kg and a cross-sectional area of {:g} m^2.'
      .format(mission.spacecraft_mass, mission.spacecraft_area))


'''
C=consumtion
N = number of particles leaving the box
n = number of engines

C = N*n*m/t, 
p = N*m*n*v
F = p/t = C*v
'''

consumption = 0.47
init_mass = 10**5
speed_boost = 1000
thrust_force = 3.13*10**3

def consumed(thrust_force,consumption,init_mass,speed_boost):
      time = speed_boost*init_mass/thrust_force
      consumed = consumption*time
      print(time/60)
      return consumed

print(consumed(thrust_force,consumption,init_mass,speed_boost))
