#EGEN KODE
#SIMPLE CONSUMPTION FUNCTION
consumption = 0.47
init_mass = 10**5
speed_boost = 1000
thrust_force = consumption*7000#consumption*velocity

def consumed(thrust_force,consumption,init_mass,speed_boost):
      time = speed_boost*init_mass/thrust_force
      consumed = consumption*time
      return consumed

print(consumed(thrust_force,consumption,init_mass,speed_boost))
