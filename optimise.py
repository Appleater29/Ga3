from hydraulic import hydraulic_cold, hydraulic_hot, hydraulic_iteration
from LMTD import log_thermal
from NTU import NTU_method, find_C_r
import numpy as np


def solution(year, Tmethod, N, N_b, passes, shell_passes, arrange, L):
    if arrange == "triangular":
        c = 0.15
        a = 0.2
    elif arrange == "sqaure":
        c = 0.2
        a = 0.34
    else:
        raise ValueError("invalid arrangement of pipes")
    mdot_1 = hydraulic_iteration(year, N, N_b, L, a, passes, shell_passes, 0.05, 0.65, side="cold")
    print("mass flow rate 1:", mdot_1)
    mdot_2 = hydraulic_iteration(year, N, N_b, L, a, passes, shell_passes, 0.05, 0.65, side="hot")
    print("mass flow rate 2:", mdot_2)
    Re_sh = hydraulic_cold(mdot_1, N, N_b, L, a, shell_passes)[1]
    Re_tube = hydraulic_hot(mdot_2, N, L, passes, shell_passes)[1]
    if Tmethod == "ntu":
        eff = NTU_method(mdot_1, mdot_2, Re_sh, Re_tube, passes, shell_passes)
    elif Tmethod == "lmtd":
        # ignore passes and shell passes for now
        eff = log_thermal(mdot_1, mdot_2, Re_sh, Re_tube, N, L, c)[4]
    else:
        raise ValueError("method not available")
    return eff
    

def find_mass(N, N_b, passes, shell_passes, L):
    L_h = 0.1 # length of headers
    m_tubes = 0.2*L*N
    m_nozzles = 0.025*4
    m_shell = 0.65*(L+L_h)
    V_endcaps = 2*np.pi*((0.07/2)**2)*0.002
    rho_plastic = 1150
    m_endcaps = V_endcaps*rho_plastic
    V_tubesheets = V_endcaps*9/2 - N*(np.pi*(0.008/2)**2)
    m_tubesheets = V_tubesheets*rho_plastic
    m_o_rings = 0.0008*2*N + 6*0.0053
    m_baffles = N_b*2.39*0.8*np.pi*(0.064/2)**2
    if passes == 1:
        m_splitters = 0
    elif passes == 2:
        m_splitters = 2.39*L_h*0.064
    elif passes == 3:
        m_splitters = 0.9*2*2.39*L_h*0.064
    elif passes > 3:
        m_splitters = 0.7*2.39*L_h*0.064*(passes-1)
    else:
        raise ValueError("invalid number of passes")
    if shell_passes > 1:
        m_splitters += L*0.064*(shell_passes-1)
    m_other = 0.01
    m = m_tubes + m_nozzles + m_shell + m_endcaps + m_tubesheets + m_o_rings + m_baffles +  m_splitters + m_other
    # print(m_baffles, m_endcaps, m_nozzles, m_o_rings, m_other, m_shell, m_splitters, m_tubes, m_tubesheets)
    # m = 0.5
    # print(m)
    return m


def test_all(L_list, N_list, N_b_list, passes_list = [1], shell_passes_list = [1], arrange = "triangular",Tmethod = "ntu", M = 1.2, year = 2025):
    design_dict = {}
    for N in N_list:
        for N_b in N_b_list:
            for passes in passes_list:
                for shell_passes in shell_passes_list:
                    for L in L_list:
                        d_parameters = (N, N_b, passes, shell_passes, arrange, L)
                        m = find_mass(N, N_b, passes, shell_passes, L)
                        if m < M:
                            outcome = solution(year, Tmethod, N, N_b, passes, shell_passes, arrange, L)
                            # find solution for effectiveness
                            design_dict.update({d_parameters: outcome})
                        elif m >= M:
                            design_dict.update({d_parameters: "too heavy"})
                        else:
                            raise ValueError("invalid mass")
    return design_dict


def find_best(design_dict):
    design_list = design_dict.items()
    print(design_list)
    best_design = 0
    for design in design_list:
        # search for highest effectiveness
        outcome = design[1]
        print(outcome)
        if outcome == "too heavy":
            print("too heavy")
        elif outcome > best_design:
            best_design = design
    return best_design


L = [264, 150, 150, 236, 236, 260, 260, 210 ,210 , 250 ,250, 163, 163]/1000
N = [12, 18, 18, 10, 10, 12, 12, 16, 16, 12, 12, 20, 20]
passes = [2, 4, 4, 2, 2, 2, 2, 1, 1, 2, 2, 4, 4]
shell_passes = [1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2]
N_b = [8, 12, 12, 8, 8, 11, 11, 6, 6, 8, 8, 6, 6]
year = [2024, 2023, 2023, 2023, 2023, 2023, 2023, 2022, 2022, 2022, 2022, 2022, 2022]
Y = [0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012]

design_dict = test_all(L, N, N_b, [4], [2], year=2023)
print(design_dict)
# print(find_best(design_dict))


# fork test


