#EGEN KODE
import numpy as np
from tqdm import trange
from numba import njit

#INITIAL VALUES
N = 10**5
L = 10**(-6)
T = 3000
n = 1000
dt = 10**(-12)
k = 1.38*10**(-23)
m = 3.3*10**(-27)

#RANDOM POSITION AND VELOCITY WITHIN BOX
def init_pos_vel(N,L,T):
    r = np.random.rand(N,3)*L
    v = np.random.normal(0, np.sqrt(k*T/m), size=(N, 3))
    return r,v

r,v = init_pos_vel(N,L,T)

#INTEGRATION LOOP
@njit
def loop(N,L,r,v,dt): 
    r = r + v*dt
    for j in range(N):
        #CHECK IF PARTICLE MOVED OUT OF BOX
        #IF SO SWITCH DIRECTION OF PARTICLE TO KEEP IN
        for k in range(3):
            if r[j,k] <= 0:
                v[j,k] = -v[j,k]

            if r[j,k] >= L:
                v[j,k] = -v[j,k]

        
    return r,v


for i in trange(n):
    r,v = loop(N,L,r,v,dt)

