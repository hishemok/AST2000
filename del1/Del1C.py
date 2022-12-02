#EGEN KODE
import numpy as np
from tqdm import trange
from numba import njit
from ast2000tools import constants

N = 10**6
L = 10**(-6)
T = 4000
n = 1000
dt = 10**(-12)
k = 1.38*10**(-23)
m = constants.m_H2#3.3*10**(-27)

def init_pos_vel(N,L,T):
    r = np.random.rand(N,3)*L
    v = np.random.normal(0, np.sqrt(k*T/m), size=(N, 3))
    return r,v

r,v = init_pos_vel(N,L,T)



@njit
def loop(N,L,r,v,dt,counter,velocity): 
    r = r + v*dt
    for j in range(N):
        #COUNT PARTICLES LEAVING THE BOX
        #STILL KEEP THEM INSIDE BUT COUNT AS IF IT LEFT
        if r[j,0] <= L/4 and r[j,1] <= L/4 and r[j,2] <= 0:
            counter +=1
            velocity += np.linalg.norm(v[j])

        
        if r[j,0] <= 0 or r[j,0] >=L:
            v[j,0] = -v[j,0]
        if r[j,1] <= 0 or r[j,1] >=L:
            v[j,1] = -v[j,1]
        if r[j,2] <= 0 or r[j,2] >=L:
            v[j,2] = -v[j,2]
        
                

        
    return r,v,counter,velocity

counter,velocity = 0,0
for i in trange(n):
    r,v,counter,velocity = loop(N,L,r,v,dt,counter,velocity)

print(counter)
print(velocity/counter)
mean_vel = velocity/counter


h2_mass = m
time = 1e-9
force = counter*mean_vel*h2_mass*1.6e13/time
consumption = counter*h2_mass*1.6e13/time
print(consumption)
print(force)
