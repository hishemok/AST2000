import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

#2.1
N = 10**5
T = 3000
k = 1.38*10**(-23)
m = 3.347*10**(-27)


def Probability(vx):
    return (m/(2*np.pi*k*T))**(1/2)*np.exp((-1/2)*((m*(vx*10**4)**2)/(k*T)))

vx =  np.linspace(-2.5, 2.5, N)

plt.subplot(2,1,1)
plt.plot(vx,Probability(vx))
plt.xlabel("V$x$")
plt.ylabel("P(v$x$)")


#2.2

def integrate(a,b):
    res,err = quad(Probability,a,b)
    return res

print(integrate(0.5,3))

#2.3
N = 10**5
T = 3000
k = 1.38*10**(-23)
m = 3.347*10**(-27)


def Probability(v):
    return (m/(2*np.pi*k*T))**(3/2)*np.exp((-1/2)*((m*(v*10**4)**2)/(k*T)))*4*np.pi*v**2

v =  np.linspace(0, 3, N)

plt.subplot(2,1,2)
plt.plot(v,Probability(v))
plt.xlabel("v")
plt.ylabel("P(v)")
plt.show()