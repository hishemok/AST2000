#EGEN KODE


import numpy as np
from ast2000tools import utils
from ast2000tools import constants
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)

angles = mission.star_direction_angles
doppler_shift = mission.star_doppler_shifts_at_sun
print(angles)

def rocket_velocity(delta_lambda_1,delta_lambda_2):
    lambda_0 = 656.3
    c = constants.c_AU_pr_yr

    v_rad_star1 = doppler_shift[0]*c/lambda_0
    v_rad_star2 = doppler_shift[1]*c/lambda_0

    v_rad_r1 = delta_lambda_1*c/lambda_0
    v_rad_r2 = delta_lambda_2*c/lambda_0

    angle1 = angles[0]*np.pi/180
    angle2 = angles[1]*np.pi/180

    k = 1/np.sin(angle2-angle1)*np.array(([np.sin(angle2),-np.sin(angle1)],[-np.cos(angle2),np.cos(angle1)]))
    x,y = np.matmul(k,np.array(([v_rad_star1-v_rad_r1],[v_rad_star2-v_rad_r2])))

    return x,y


dop_shift = [-0.0353227073295147, -0.0277205765307337]
#Velocity: (-2.15192, 2.00343) AU/yr

vx,vy = rocket_velocity(dop_shift[0],dop_shift[1])
print(vx,vy)
#[-2.15192058] [2.00342743]