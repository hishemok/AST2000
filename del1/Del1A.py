#EGEN KODE
from scipy.integrate import quad
import numpy as np

sigma = 1
my = 0
def f(x):
    return 1/(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5*((x-my)/sigma)**2)


def integrate(a,b):
    res,err = quad(f,a,b)
    return res

print(integrate(-sigma,sigma))
print(integrate(-2*sigma,2*sigma))
print(integrate(-3*sigma,3*sigma))

import matplotlib.pyplot as plt
x = np.linspace(-3*sigma,3*sigma,1000)
plt.plot(x,f(x))
plt.hlines(f(0)/2,xmin=-3*sigma,xmax=3*sigma)
plt.show()

