import numpy as np
from tqdm import trange
from numba import njit

N = 10**5
L = 10**(-6)
T = 3000
n = 1000
dt = 10**(-12)
k = 1.38*10**(-23)
m = 3.3*10**(-27)

def init_pos_vel(N,L,T):
    r = np.random.rand(N,3)*L
    v = np.random.normal(0, np.sqrt(k*T/m), size=(N, 3))
    return r,v

r,v = init_pos_vel(N,L,T)

@njit
def loop(N,L,r,v,dt,counter,velocity): 
    r = r + v*dt
    for j in range(N):
        if r[j,0] <= L/4 and r[j,1] <= L/4 and r[j,2] <= 0:
            counter +=1
            velocity += np.linalg.norm(v[j])
        for k in range(3):
            if r[j,k] <= 0:
                v[j,k] = -v[j,k]

            if r[j,k] >= L:
                v[j,k] = -v[j,k]

        
    return r,v,counter,velocity

counter,velocity = 0,0
for i in trange(n):
    r,v,counter,velocity = loop(N,L,r,v,dt,counter,velocity)

print(counter)
print(velocity/counter)
mean_vel = velocity/counter
print(f'Force:{counter*mean_vel*m*1.6*10**13*10**9}')
