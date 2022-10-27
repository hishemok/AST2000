#EGEN KODE

import matplotlib.pyplot as plt
import numpy as np
from ast2000tools import utils
from ast2000tools import constants
utils.check_for_newer_version()
seed = utils.get_seed('<Claudieg>')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)


planet_positions0 = np.einsum('ij->ji',system.initial_positions)
planet_trajectories = np.load('planet_trajectories.npz')
planet_positions = np.einsum('ijk->kji',planet_trajectories['planet_positions'])[:,2:5]

#mission.measure distances 
distance_list = [1.29816269e-04, 1.33599513e+01, 9.30388861e+00, 3.84511610e+01,1.08823373e+01, 2.19127986e+01, 6.45994960e+00, 6.66324004e+00]






def find_rocket(time,distances):
    positions = planet_positions[time]
    #assume distance contains information about 3 planets
    r1 = distances[0]
    r2 = distances[1]
    r3 = distances[2]
    x1,y1 = positions[0,0],positions[0,1]
    x2,y2 = positions[1,0],positions[1,1]
    x3,y3 = positions[2,0],positions[2,1]

    
    a = -2*x1 + 2*x2
    b = -2*y1 + 2*y2
    c = r1**2 - r2**2 - x1**2 + x2**2 - y1**2 + y2**2
    d = -2*x2 + 2*x3
    e = -2*y2 + 2*y3
    f = r2**2 - r3**2 - x2**2 + x3**2 - y2**2 + y3**2

    x = (c*e-f*b)/(e*a-b*d)
    y = (c*d-a*f)/(b*d-a*e)

    #test
    
    angle = np.linspace(0,2*np.pi,1000)
    plt.plot(x1+np.cos(angle)*r1,y1+np.sin(angle)*r1,label='1')
    plt.plot(x2+np.cos(angle)*r2,y2+np.sin(angle)*r2,label='2')
    plt.plot(x3+np.cos(angle)*r3,y3+np.sin(angle)*r3,label='3')
    plt.scatter(x,y)
    plt.legend()
    plt.axis('equal')
    plt.show()
    
    
    return np.array([x,y])

position = find_rocket(15993,distance_list[2:5])
print(position)



