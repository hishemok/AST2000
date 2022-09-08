import numpy as np

consumption = 1500 #consumption/0.33 (where 0.33 is C by one sheet) engines stacked
init_mass = 5*10**5
speed_boost = 1000
thrust_force = consumption*7800 #init_mass*30#
print(thrust_force)

def consumed(thrust_force,consumption,init_mass,speed_boost):
      time = speed_boost*init_mass/thrust_force
      consumed = consumption*time
      time = time/60
      return consumed,time




r0 =  8552.29*10**3

v_esc = 14270.788097398781 #np.sqrt(2*g*planet_mass/r0)
def gravity(r,init_mass):
    planet_mass = 1.3056371666590186e+25
    g = 6.67*10**(-11)
    return init_mass*g*planet_mass/r**2


total_time = 0
i = 0
while i*speed_boost <= v_esc:
    Force = thrust_force-gravity(r0,init_mass)
    C,time = consumed(Force,consumption,init_mass,speed_boost)
    r0 = r0 + speed_boost
    init_mass -= C
    total_time += time
    i+=1
print(f'Remaining mass:{init_mass:.2f}kg, time:{total_time:.2f}min')
print(i*speed_boost)



