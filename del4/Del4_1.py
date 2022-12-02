#EGEN KODE

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from tqdm import trange
from ast2000tools import utils
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)

img = Image.open('AST2000/sample0000.png')
pixels = np.array(img) # png into numpy array
width = len(pixels[0, :])
height = len(pixels[:,0])
print(width,height)
#640x480


alpha = 70*np.pi/180
theta_0 = np.pi/2
phi_0 = 0
def create_grid(alpha,theta_0,phi_0):
    #follows given method to create a grid to set rbg colors in.
    xmin = -2*np.sin(alpha/2)/(1+np.cos(alpha/2))
    xmax = 2*np.sin(alpha/2)/(1+np.cos(alpha/2))
    ymin = -2*np.sin(alpha/2)/(1+np.cos(alpha/2))
    ymax = 2*np.sin(alpha/2)/(1+np.cos(alpha/2))

    x = np.linspace(xmin,xmax,width)
    y = np.linspace(ymin,ymax,height)

    X,Y = np.meshgrid(x,y)

    rho = np.sqrt(X**2+Y**2)
    beta = 2*np.arctan(rho/2)

    
    theta = theta_0 + np.arcsin(np.cos(beta)*np.cos(theta_0) + Y/rho*np.sin(beta)*np.sin(theta_0))

    
    phi = phi_0 + np.arctan((X*np.sin(beta)/(rho*np.sin(theta_0)*np.cos(beta)-Y*np.cos(theta_0)*np.sin(beta))))

    
    grid = np.zeros((height,width,3),dtype = "uint8")
    return grid, theta, phi


grid,theta,phi = create_grid(alpha,theta_0,phi_0)


himmelkule = np.load('Downloads/himmelkule.npy')

#uses angle to draw right rgb from "himmelkule" to right angles
for i in range(height):
    for j in range(width):
        pixel_index = mission.get_sky_image_pixel(theta[i,j],phi[i,j])
        dont,need,r,g,b = himmelkule[pixel_index]
        grid[i,j] = r,g,b

img2 = Image.fromarray(grid)
img2.save('AST2000/picture.png') 

#Does the same but from angles 1-360 for phi and saves pictures individually.
for angle in trange(1,361):
    alpha = 70*np.pi/180
    theta_0 = np.pi/2
    phi_0 = angle*np.pi/180
    grid,theta,phi = create_grid(alpha,theta_0,phi_0)
    for i in range(height):
        for j in range(width):

            pixel_index = mission.get_sky_image_pixel(theta[i,j],phi[i,j])
            dont,need,r,g,b = himmelkule[pixel_index]
            grid[i,j] = r,g,b


    img2 = Image.fromarray(grid)
    img2.save(f'AST2000/sky_image_angle_{angle}.png') 

