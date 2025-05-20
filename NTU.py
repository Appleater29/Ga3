import numpy as np
from LMTD import find_H
from fixed_constants import *


# N = 13
# L = 0.35

def find_C_r(mdot_1, mdot_2):
    C_min = min(mdot_1, mdot_2) * cp
    C_max = max(mdot_1, mdot_2) * cp
    return C_min / C_max

def NTU_method(mdot_1, mdot_2, Re_sh, Re_tube, N, L,c = 0.15, passes = 1, shell_passes = 1):
    C_min = min(mdot_1, mdot_2) * cp
    C_max = max(mdot_1, mdot_2) * cp
    C_r = C_min / C_max
    H = find_H(Re_sh, Re_tube, L, c)
    A = (np.pi)*N*d_i*L
    ntu = H*A/C_min
    if C_r == 1 and passes == 1:
        eff = ntu/(1+ntu) # for balanced counter flow
    elif C_r != 1 and passes == 1:
        eff = (1-np.exp(-ntu*(1-C_r)))/(1-C_r*np.exp(-ntu*(1-C_r))) # for normal counter flow
    elif passes%2 == 0 and shell_passes == 1:
        eff = 2/(1+C_r+((1+C_r)**0.5)*(1+np.exp(-ntu*(1+C_r**2)**0.5))/(1-np.exp(-ntu*(1+C_r**2)**0.5)))
    elif passes%2 == 0:
        eff1 = 2/(1+C_r+((1+C_r)**0.5)*(1+np.exp(-ntu*(1+C_r**2)**0.5))/(1-np.exp(-ntu*(1+C_r**2)**0.5)))
        eff = (((1-eff1*C_r)/(1-eff1))**shell_passes - 1)/(((1-eff1*C_r)/(1-eff1))**shell_passes -C_r)
    else:
        raise ValueError("Incorrect NTU parameters")
    return eff

# print(NTU_method(0.5, 0.47, 15281, 11283*0.47/0.45, N, L))