import numpy as np
# fixed variables
# T1_in = 20 # degrees
T1_in = 19.8 # degrees
# T2_in = 60 # degrees 
T2_in = 51.7 # degrees 
cp = 4179 # J/kgK
rho_w = 990.1 # kg/ms
k_w = 0.632 # W/mK
mu = 6.51e-4 # kg/ms
Pr = 4.31 
k_tube = 386 # W/mK
d_o = 0.008 # meters
d_i = 0.006 # meters
d_noz = 0.02 # meters
d_sh = 0.064 # meters
# L = 0.35 # meters
Y = 0.014 # meters

# constant areas
A_tube = np.pi * 0.25 * d_i**2
A_noz = np.pi * 0.25 * d_noz**2
A_pipe = np.pi * 0.25 * d_sh**2
A_hose = np.pi * 0.25 * 0.025**2
