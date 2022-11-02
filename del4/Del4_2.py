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


def create_grid(alpha,theta_0,phi_0,height=480,width=640):
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


def find_angle_of_image(img0,img2):
    himmelkule = np.load('/Users/hishem/AST2000/himmelkule.npy')

    pixels = np.array(img0)
    width = len(pixels[0, :])
    height = len(pixels[:,0])

    alpha = 70*np.pi/180
    theta_0 = np.pi/2
    max_angle = 361
    #creates photos for all phi between 1 and 360 to compare with image at phi = 0
    #If the angle between the created image and the original is equal to the angle between the given image
    # and the original, then the loop stops and the function returns what angle the photo is at.
    for angle in trange(1,max_angle):
        
        phi_0 = angle*np.pi/180
        grid,theta,phi = create_grid(alpha,theta_0,phi_0)

        err = np.zeros(max_angle)
        img2_pixels = np.array(img2)
        corr = np.sum(pixels-img2_pixels)
        for i in range(height):
            for j in range(width):
                pixel_index = mission.get_sky_image_pixel(theta[i,j],phi[i,j])
                dont,need,r,g,b = himmelkule[pixel_index]
                grid[i,j] = r,g,b

                if i == (height-1) and j == (width-1):
                    #when the picture at phi = angle is created i compare the angles and checks if its the right angle.
                    #if yes then we break the loop. if no then we move to the next iteration.
                    err[angle] = np.sum(pixels - grid)
                    if err[angle] == corr:
                        return angle

#sample image and a random image is given to find the angle between them.
img0 = Image.open('/Users/hishem/AST2000/del4/sample0000.png')
img2 = Image.open('/Users/hishem/AST2000/sky_picture.png')


print('Angle of image:',find_angle_of_image(img0,img2,))
#Angle of image: 199

