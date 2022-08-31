import numpy as np


L = 10**(-6)
T = 3*10**3
N = 100
n = 100
dt = 0.1

def init(N,L,T):
    r = np.random.rand((N,3))*L
    v = np.random.rand((N,3))*L
    return r,v

r,v = init(N,L,T)

for i in range(n):
    r += 
