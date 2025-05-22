import numpy as np
from fixed_constants import *

def hydraulic_cold_old(m1_dot, N, N_b, L, shell_passes, layout="tri"):
    if layout == 'tri':
        a = 0.2
    elif layout == 'square':
        a = 0.34
    
    B = L / (N_b + 1)
    A_sh = d_sh * (Y - d_o) * B / (Y*shell_passes)
    # cold stream
    v_sh = m1_dot / (rho_w * A_sh)
    d_sh_dash = d_sh * A_sh / A_pipe 
    Re_sh = rho_w * v_sh * d_sh_dash / mu
    delta_p_sh = 4 * a * np.power(Re_sh, -0.15) * N * rho_w * v_sh**2

    v_noz1 = m1_dot / (rho_w * A_noz)
    delta_p_noz1 = rho_w * v_noz1**2

    v_hose_max = 0.658 / (rho_w * A_hose)
    k_hose = 15840 / (0.5 * rho_w * v_hose_max**2)
    v_hose = m1_dot / (rho_w * A_hose)
    delta_p_hose = 0.5 * rho_w * v_hose**2 * k_hose

    delta_p_1 = delta_p_sh + delta_p_noz1 + delta_p_hose
    print(delta_p_1)

    return(delta_p_1, Re_sh)

def hydraulic_cold_Kern(m1_dot, N, N_b, L, shell_passes, layout="tri"):
    # cold stream
    B = L / (N_b + 1)
    A_sh = d_sh * (Y - d_o) * B / (Y*shell_passes) # not sure if shell passes are needed
    v_sh = m1_dot / (rho_w * A_sh)

    if layout == "square":
        D_e = (4 * (Y**2 - (np.pi * d_o**2 / 4))) / (np.pi * d_o)
    elif layout == "tri":
        D_e = (4 * ((Y**2 * np.sqrt(3) / 4) - (np.pi * d_o**2 / 8))) / (np.pi * d_o / 2)
    else:
        raise ValueError("Layout must be 'square' or 'triangular'")
    
    Re_sh = rho_w * v_sh * D_e / mu

    G_s = m1_dot / A_sh
    f = np.exp(0.576 - 0.19 * np.log(Re_sh))
    delta_p_sh = f * G_s**2 * (N_b  + 1) * d_sh / (2 * rho_w * D_e)

    v_noz1 = m1_dot / (rho_w * A_noz)
    delta_p_noz1 = rho_w * v_noz1**2

    v_hose_max = 0.658 / (rho_w * A_hose)
    k_hose = 15840 / (0.5 * rho_w * v_hose_max**2)
    v_hose = m1_dot / (rho_w * A_hose)
    delta_p_hose = 0.5 * rho_w * v_hose**2 * k_hose
    
    delta_p_1 = delta_p_sh + delta_p_noz1 + delta_p_hose
    print(delta_p_1)

    return(delta_p_1, Re_sh)

def hydraulic_hot_old(m2_dot, N, L, passes, shell_passes):
    # hot stream
    sigma = N * (A_tube/passes) / (A_pipe/shell_passes)
    mtube_dot = m2_dot / N
    v_tube = mtube_dot / (rho_w * A_tube)
    Re_tube = rho_w * v_tube * d_i / mu

    v_noz2 = m2_dot / (rho_w * A_noz)

    delta_p_tube = 0.5 * rho_w * v_tube**2 * friction_factor(Re_tube) * L / d_i
    delta_p_ends = 0.5 * rho_w * v_tube**2 * (K_c(sigma) + K_e(sigma))
    delta_p_noz2 = rho_w * v_noz2**2

    v_hose_max = 0.436 / (rho_w * A_hose)
    k_hose = 9320 / (0.5 * rho_w * v_hose_max**2)
    v_hose = m2_dot / (rho_w * A_hose)
    delta_p_hose = 0.5 * rho_w * v_hose**2 * k_hose

    delta_p_2 = (delta_p_tube + delta_p_ends) * passes + delta_p_noz2 + delta_p_hose
    print(delta_p_2)

    return(delta_p_2, Re_tube)

def hydraulic_hot_Kern(m2_dot, N, L, passes, shell_passes):
    # hot stream
    mtube_dot = m2_dot / N
    v_tube = mtube_dot / (rho_w * A_tube)
    Re_tube = rho_w * v_tube * d_i / mu

    v_noz2 = m2_dot / (rho_w * A_noz)

    f = (1.58*np.log(Re_tube) - 3.28)**(-2)
    delta_p_tot = 2*passes*(f*L/d_i + 1)*(rho_w*v_tube**2)
    delta_p_noz2 = rho_w * v_noz2**2

    v_hose_max = 0.436 / (rho_w * A_hose)
    k_hose = 9320 / (0.5 * rho_w * v_hose_max**2)
    v_hose = m2_dot / (rho_w * A_hose)
    delta_p_hose = 0.5 * rho_w * v_hose**2 * k_hose

    delta_p_2 = delta_p_tot + delta_p_noz2 + delta_p_hose
    print(delta_p_2)

    return(delta_p_2, Re_tube)

def friction_factor(Re):
    return(np.power((1.82*np.log10(Re) - 1.64), -2))

def K_c(x):
    return(0.5 - 0.4*x)

def K_e(x):
    if(x < 0.5):
        return(x**2 - 2.1 * x + 1)
    else:
        raise ValueError("sigma too high, check number of tubes")


def cold_chic(x, year=2025):
    """Return compressor cold-side pressure drop (Pa) for given mass flow rate x (L/s), based on year."""
    if year == 2022:
        return 1e5 * (-0.7843 * x**2 - 0.4802 * x + 0.6598)
    elif year == 2023:
        return 1e5 * (-0.4914 * x**2 - 0.3966 * x + 0.6612)
    elif year == 2024:
        return 1e5 * (-0.7691 * x**2 - 0.3763 * x + 0.6555)
    elif year == 2025:
        return 1e5 * (-0.5557 * x**2 - 0.4559 * x + 0.7049)
    else:
        raise ValueError("Unsupported year for cold_chic. Choose 2022, 2023, 2024, or 2025.")


def hot_chic(x, year=2025):
    """Return compressor hot-side pressure drop (Pa) for given mass flow rate x (L/s), based on year."""
    if year == 2022:
        return 1e5 * (-0.4713 * x**2 - 0.8896 * x + 0.6381)
    elif year == 2023:
        return 1e5 * (-0.8204 * x**2 - 0.6777 * x + 0.5626)
    elif year == 2024:
        return 1e5 * (-0.8117 * x**2 - 0.6792 * x + 0.6145)
    elif year == 2025:
        return 1e5 * (-1.2312 * x**2 - 0.6973 * x + 0.6251)
    else:
        raise ValueError("Unsupported year for hot_chic. Choose 2022, 2023, 2024, or 2025.")
    
def hydraulic_iteration(year, N, N_b, L, a, passes, shell_passes, x_min, x_max, tol=0.001, side="cold"):
    
    B = L / (N_b + 1)

    # Define the appropriate error function
    if side == "cold":
        def error(x):
            delta_p = hydraulic_cold_Kern(x, N, N_b, L, a, shell_passes)[0]
            return(delta_p - cold_chic(x, year))
    elif side == "hot":
        def error(x):
            delta_p = hydraulic_hot_Kern(x, N, L, passes, shell_passes)[0]
            return(delta_p - hot_chic(x, year))
    else:
        raise ValueError("Invalid side. Use 'cold' or 'hot'.")

    # Check if initial bounds bracket a root
    f_min = error(x_min)
    f_max = error(x_max)

    if f_min * f_max > 0:
        raise ValueError("No root in interval: error has same sign at both ends")

    # Bisection loop
    while (x_max - x_min) > tol:
        x_mid = 0.5 * (x_min + x_max)
        f_mid = error(x_mid)

        if f_mid * f_min > 0 and f_mid * f_max > 0:
            raise ValueError("Two roots or no sign change in interval")

        elif f_mid * f_min > 0:
            x_min = x_mid
            f_min = f_mid
        elif f_mid * f_max > 0:
            x_max = x_mid
            f_max = f_mid
        else:
            return x_mid  # Exact root found
    delta_p = hydraulic_cold_Kern(x_mid, N, N_b, L, a, shell_passes)[0]
    print("hot pressure:", delta_p)

    return 0.5 * (x_min + x_max)

# m1_dot = hydraulic_iteration(0.05, 0.65, side="cold")
# m2_dot = hydraulic_iteration(0.05, 0.65, side="hot")

# delta_p_1, Re_sh = hydraulic_cold(m1_dot)
# delta_p_2, Re_tube = hydraulic_hot(m2_dot)

# print(m1_dot, m2_dot, Re_sh, Re_tube, delta_p_1, delta_p_2)