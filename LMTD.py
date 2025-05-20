import numpy as np
import matplotlib.pyplot as plt

from fixed_constants import *

# N = 13
# N_b = 9
# d_o = 0.008 # meters
# d_i = 0.006 # meters
# d_noz = 0.02 # meters
# d_sh = 0.064 # meters
# L = 0.35 # meters
# Y = 0.014 # meters
# T1_in = 20 #degrees
# T2_in = 60 #degrees 
# cp = 4179 #J/kgK
# rho_w = 990.1 #kg/ms
# k_w = 0.632 #W/mK
# mu = 6.51e-4 #kg/ms
# Pr = 4.31 
# k_tube = 386 #W/mK
# F = 1
# c = 0.15 # square to match example

def find_H (Re_sh, Re_tube, L, c = 0.15 ):
    Nu_i = 0.023*((Re_tube)**0.8)*(Pr**0.3)
    Nu_o = c*(Re_sh**0.6)*(Pr**0.3)
    h_i = Nu_i*k_w/d_i
    h_o = Nu_o*k_w/d_o
    A_i = (np.pi)*d_i*L
    A_o = (np.pi)*d_o*L
    return (1/h_i +A_i*(np.log(d_o/d_i))/(2*(np.pi)*k_tube*L) + (1/h_o)*(A_i/A_o))**(-1)

def logmtemp(T1_out, T2_out):
    T_lm = ((T2_in - T1_out) - (T2_out - T1_in)) / (np.log((T2_in - T1_out)/(T2_out - T1_in)))
    return T_lm


def correction_factor(T1_out, shell_passes):
    x = (T1_out - T1_in) / (T2_in - T1_in)

    # assume R = 1.5
    if shell_passes == 1:
        p_vals = [0, 0.2, 0.3, 0.35, 0.4, 0.45, 0.46]
        f_vals = [1, 0.98, 0.95, 0.92, 0.8, 0.6, 0.0]
    elif shell_passes == 2:
        p_vals = [0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.57]
        f_vals = [1, 0.99, 0.98, 0.96, 0.93, 0.85, 0.7, 0.0]
    else:
        raise ValueError("Only 1 or 2 shell passes supported")

    # Handle out-of-range
    if x <= p_vals[0]:
        return f_vals[0]
    if x >= p_vals[-1]:
        return f_vals[-1]

    # Linear interpolation
    for i in range(1, len(p_vals)):
        if x < p_vals[i]:
            x0, x1 = p_vals[i - 1], p_vals[i]
            y0, y1 = f_vals[i - 1], f_vals[i]
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0) 

def fT1_out(T1_out, H, A, mdot_1, mdot_2):
        F = 1
        CF = correction_factor(T1_out)
        T2_out = T2_in - (mdot_1/mdot_2)*(T1_out-T1_in)
        T_lm = logmtemp(T1_out, T2_out)
        return mdot_1*cp*(T1_out-T1_in)-H*A*T_lm*F


def log_thermal(mdot_1, mdot_2, Re_sh, Re_tube, N, L, c=0.15):
    A = (np.pi)*N*d_i*L
    H = find_H(Re_sh, Re_tube, L, c)
    T1_out_max = T2_in - 5
    # T1_out_max = T2_in - 1
    T1_out_min = T1_in + 5
    # T1_out_min = T1_in + 1
    
    # find T1_out by iteration:

    if T1_out_max == T1_out_min:
        print("Temperatures too close")

    while T1_out_max - T1_out_min > 0.0001:
        T1_out_mid = (T1_out_min + T1_out_max)/2  
        f_max = fT1_out(T1_out_max, H, A, mdot_1, mdot_2)
        f_min = fT1_out(T1_out_min, H, A, mdot_1, mdot_2)
        f_mid = fT1_out(T1_out_mid, H, A, mdot_1, mdot_2)
        if np.sign(f_max) == np.sign(f_min):
             print("error - no root in interval")
             break
        elif np.sign(f_mid) == np.sign(f_min):
             T1_out_min = T1_out_mid
        elif np.sign(f_max) == np.sign(f_mid):
             T1_out_max = T1_out_mid
        else:
             print("error - two roots in interval")
             break

    T1_out = T1_out_mid
    T2_out = T2_in - (mdot_1/mdot_2)*(T1_out-T1_in)
    T_lm = logmtemp(T1_out, T2_out)
    Qdot = mdot_1*cp*(T1_out - T1_in)
    eff = Qdot/(mdot_2 * cp*(T2_in - T1_in))

    # plotting
    # x = np.linspace(-100, 100, 1000)
    # y_2 = np.zeros(1000)
    # y = fT1_out(x, H, (np.pi)*N*d_i*L, mdot_1, mdot_2)
    # plt.plot(x,y)
    # plt.plot(x,y_2)
    # plt.show()

    return T1_out, T2_out, T_lm, Qdot, eff

# print(log_thermal(0.5, 0.47, 15281, 11283*0.47/0.45))

# x = np.linspace(-100, 100, 1000)
# y_2 = np.zeros(1000)
# y = fT1_out(x, 3600, (np.pi)*N*d_i*L, 0.50, 0.47)
# plt.plot(x,y)
# plt.plot(x,y_2)
# plt.show()