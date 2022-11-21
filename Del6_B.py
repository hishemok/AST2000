#EGEN KODE
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from numba import njit
from ast2000tools import utils
from ast2000tools import constants
utils.check_for_newer_version()
seed = utils.get_seed('Claudieg')
from ast2000tools import space_mission 
mission = space_mission.SpaceMission(seed)
from ast2000tools import solar_system
system = solar_system.SolarSystem(seed)
c = constants.c
k = constants.k_B

'''Max doppler shift'''
def doppler_shift(lambda0):
    return lambda0*10*1e3/c
''''''
'''expression for sigma'''
def sigma(lambda0,T,m):
    return (lambda0/c)*np.sqrt(k*T/m)
''''''

spectrum = np.load('spectrum.npy')
sigma_noise = np.load('sigma_noise.npy')

def F(Fmin,lam,lam0,T,m):
    return 1 + (Fmin-1)*np.exp(-0.5*((lam-lam0)/sigma(lam0,T,m))**2)


def Spectral_Analysis(spectrum,sigma_noise,gass):
    u = 1.066*10**(-27)
    m = float(gass[0])*u

    lam0 = int(gass[1])
    delta_lambda = doppler_shift(lam0)
    tol = 1e-3
    
    #use doppler shift delta_lambda to know upper and lower limit of where to look
    lower_idx = np.where(abs(spectrum[:,0] - np.round(lam0 - delta_lambda, 3)) < tol)[0][0]
    upper_idx = np.where(abs(spectrum[:,0] - np.round(lam0 + delta_lambda, 4)) < tol)[0][-1]

    #cut spectrum with new lower and upper indexes
    #also *=1e-9 to calculate with [m]
    spectrum = spectrum[lower_idx:upper_idx]
    spectrum[:,0] *= 1e-9
    sigma_noise = sigma_noise[lower_idx:upper_idx]
    sigma_noise[:,0] *= 1e-9
    
    T = np.linspace(150, 450, 100)
    Fmin = np.linspace(0.7, 1, 10)
    Lambda = np.linspace(spectrum[0,0],spectrum[-1,0],50)
    #big number, but will most likely get new value in first loop. 
    min = 1e9
    for i in trange(len(T)):
        for f in Fmin:
            for l in Lambda:
                calc = F(f,l,spectrum[:,0],T[i],m)
                chi_sq = np.sum( ( (spectrum[:,1]-calc) /sigma_noise[:,1])**2 )

                if chi_sq < min and f != 0.7 and f != 1:
                    min = chi_sq
                    #save for which values gave the smallest chi_sq
                    vals = np.array([f,l*1e9,T[i]])

    spectrum[:,0] = spectrum[:,0]*1e9 #back to original size
    plt.plot(spectrum[:,0],spectrum[:,1],label='Measured data')
    plt.xlabel('wavelength nm')
    plt.ylabel('flux')
    plt.title(gass[2])
    plt.plot(spectrum[:,0],F(vals[0],vals[1],spectrum[:,0],vals[2],m),label='Gaussian line profile')
    plt.legend()
    plt.show()
    return vals



all_values = np.zeros((12,5))
O2_632 = np.array([31.999,632,'O2_632'])
vals = Spectral_Analysis(spectrum,sigma_noise,O2_632)
all_values[0] = np.array([632,vals[0],vals[1],vals[2],abs(632-vals[1])])

O2_690 = np.array([31.999,690,'O2_690'])
vals = Spectral_Analysis(spectrum,sigma_noise,O2_690)
all_values[1] = np.array([690,vals[0],vals[1],vals[2],abs(690-vals[1])])

O2_760 = np.array([31.999,760,'O2_760'])
vals = Spectral_Analysis(spectrum,sigma_noise,O2_760)
all_values[2] = np.array([760,vals[0],vals[1],vals[2],abs(760-vals[1])])

H2O_720 = np.array([18.01528,720,'H2O_720'])
vals = Spectral_Analysis(spectrum,sigma_noise,H2O_720)
all_values[3] = np.array([720,vals[0],vals[1],vals[2],abs(720-vals[1])])

H2O_820 = np.array([18.01528,820,'H2O_820'])
vals = Spectral_Analysis(spectrum,sigma_noise,H2O_820)
all_values[4] = np.array([820,vals[0],vals[1],vals[2],abs(820-vals[1])])

H2O_940 = np.array([18.01528,940,'H2O_940'])
vals =Spectral_Analysis(spectrum,sigma_noise,H2O_940)
all_values[5] = np.array([940,vals[0],vals[1],vals[2],abs(940-vals[1])])

CO2_1400 = np.array([44.01,1400,'CO2_1400'])
vals = Spectral_Analysis(spectrum,sigma_noise,CO2_1400)
all_values[6] = np.array([1400,vals[0],vals[1],vals[2],abs(1400-vals[1])])

CO2_1600 = np.array([44.01,1600,'CO2_1600'])
vals = Spectral_Analysis(spectrum,sigma_noise,CO2_1600)
all_values[7] = np.array([1600,vals[0],vals[1],vals[2],abs(1600-vals[1])])

CH4_1660 = np.array([16.04,1660,'CH4_1660'])
vals = Spectral_Analysis(spectrum,sigma_noise,CH4_1660)
all_values[8] = np.array([1660,vals[0],vals[1],vals[2],abs(1660-vals[1])])

CH4_2200 = np.array([16.04,2200,'CH4_2200'])
vals = Spectral_Analysis(spectrum,sigma_noise,CH4_2200)
all_values[9] = np.array([2200,vals[0],vals[1],vals[2],abs(2200-vals[1])])

CO = np.array([28.01,2340,'CO_2340'])
vals = Spectral_Analysis(spectrum,sigma_noise,CO)
all_values[10] = np.array([2340,vals[0],vals[1],vals[2],abs(2340-vals[1])])

N2O = np.array([44.013,2870,'N2O_2870'])
vals = Spectral_Analysis(spectrum,sigma_noise,N2O)
all_values[11] = np.array([2870,vals[0],vals[1],vals[2],abs(2870-vals[1])])


print('Wavelength | Fmin | Temperature | Delta Lambda')

for i in range(11):
    print(f'{all_values[i,0]}       {all_values[i,1]:.3f}      {all_values[i,3]:.2f}       {all_values[i,4]:.4}')

'''
Wavelength | Fmin | Temperature | Delta Lambda
632.0       0.933      150.00       0.007657
690.0       0.900      340.91       0.01425
760.0       0.733      150.00       0.02283
720.0       0.967      150.00       0.01685
820.0       0.833      150.00       0.02336
940.0       0.867      150.00       0.02809
1400.0       0.833      150.00       0.04003
1600.0       0.867      237.88       0.04749
1660.0       0.933      450.00       0.05372
2200.0       0.800      168.18       0.06495
2340.0       0.767      168.18       0.06906
'''
