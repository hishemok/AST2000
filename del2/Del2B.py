import numpy as np
import matplotlib.pyplot as plt


from ast2000tools import utils
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)


n = 10**4
v = np.zeros((n,2))
r = np.zeros_like(v)

v[0] = np.einsum('ij->ji',system.initial_velocities)[0]
r[0] = np.einsum('ij->ji',system.initial_positions)[0]

dt = 1e-5
star_mass = system.star_mass
G = 4*np.pi**2
#planet 0 position integration loop
for i in range(n-1): 
    r_norm = np.linalg.norm(r[i])

    a = -G*star_mass*(r[i]/r_norm**3)
    vi05 = v[i] + a*dt/2
    r[i+1] = r[i] + vi05*dt

    r_norm = np.linalg.norm(r[i+1])

    a =  -G*star_mass*(r[i+1]/r_norm**3)
    v[i+1] = vi05 + a*dt




def A(v1,v2):#function calculates area swept after planet moved
    #calculates the angles between to position vectors
    unit_v1 = v1/np.linalg.norm(v1)
    unit_v2 = v2/np.linalg.norm(v2)
    dot_product = np.dot(unit_v1, unit_v2)
    angle = np.arccos(dot_product)
    A = 0.5*np.linalg.norm(v1)**2*angle
    ds = np.linalg.norm(v1)*angle
    return f'Area: {A}, Angle: {angle}, ds: {ds}'#A,angle,ds

tol = 0.5e-3
ry = r[:,1]
rx = r[:,0]
DT = 5

#find area swept when rx is greates +- 5 timesteps
for i in range(n-1):
    if rx[i] == np.max(rx):
        break
r1 = np.array([r[i-DT,0],r[i-DT,1]])
r2 = np.array([r[i+DT,0],r[i+DT,1]])

#scatter points to check if calculations look right
plt.scatter(r1[0],r1[1])
plt.scatter(r2[0],r2[1])
plt.hlines(r[i,1],0,r[i,0])
print(A(r1,r2))

#mean vel between the two points
mean_vel_a = np.sum(v[i:i+DT])/DT


#find area swept when ry is greates +- 5 timesteps
for i in range(n-1):
    if ry[i] == np.max(ry):
        break
r1 = np.array([r[i-DT,0],r[i-DT,1]])
r2 = np.array([r[i+DT,0],r[i+DT,1]])

#scatter points to check if calculations look right
plt.scatter(r1[0],r1[1])
plt.scatter(r2[0],r2[1])
plt.vlines(r[i,0],0,r[i,1])
print(A(r1,r2))

#mean vel between the two points
mean_vel_b = np.sum(v[i:i+DT])/DT

print(mean_vel_a)
print(mean_vel_b)


plt.plot(rx,ry)
plt.show()


# mean_vel_total = np.mean(v[:])
# print(mean_vel_total)
# abs_mean_vel_total = np.mean(np.abs(v[:]))
# print(abs_mean_vel_total)


