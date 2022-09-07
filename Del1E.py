import numpy as np

consumption = 100 #assumtion
init_mass = 10**5
speed_boost = 1000
thrust_force = 3*10**6#Estimated thrust force must be greater than gravitational force

def consumed(thrust_force,consumption,init_mass,speed_boost):
      time = speed_boost*init_mass/thrust_force
      consumed = consumption*time
      time = time/60
      return consumed,time


v_esc = 15*10**3 #approx estimated compared to earth

r0 =  9364.77
def gravity(r,init_mass):
    return np.sqrt(6.67*10**(-11)*init_mass*5.6*10**24/r**2)


total_time = 0
i = 0
while i*speed_boost <= v_esc:
    Force = thrust_force-gravity(r0,init_mass)
    C,time = consumed(Force,consumption,init_mass,speed_boost)
    r0 = r0 + speed_boost
    init_mass -= C
    total_time += time
    i+=1
print(f'Mass:{init_mass:.2f}kg, time:{total_time:.2f}min')
