from hydraulic import hydraulic_cold_Kern, hydraulic_hot_Kern, hydraulic_iteration
from LMTD import log_thermal
from NTU import NTU_method
import numpy as np
import matplotlib.pyplot as plt


def solution(year, Tmethod, N, N_b, passes, shell_passes, arrange, L):
    if arrange == "tri":
        c = 0.15
        a = 0.2
    elif arrange == "sqaure":
        c = 0.2
        a = 0.34
    else:
        raise ValueError("invalid arrangement of pipes")
    mdot_1 = hydraulic_iteration(year, N, N_b, L, a, passes, shell_passes, 0.05, 0.65, side="cold")
    mdot_2 = hydraulic_iteration(year, N, N_b, L, a, passes, shell_passes, 0.05, 0.65, side="hot")
    Re_sh = hydraulic_cold_Kern(mdot_1, N, N_b, L, shell_passes, layout=arrange)[1]
    Re_tube = hydraulic_hot_Kern(mdot_2, N, L, passes, shell_passes)[1]
    if Tmethod == "ntu":
        outcome = NTU_method(mdot_1, mdot_2, Re_sh, Re_tube, N, L, c, passes, shell_passes)
        eff = outcome[0]
        Qdot = outcome[1]

    elif Tmethod == "lmtd":
        outcome = log_thermal(mdot_1, mdot_2, Re_sh, Re_tube, N, L, c, shell_passes)
        eff = outcome[4]
        Qdot = outcome[3]
    else:
        raise ValueError("method not available")
    return eff, Qdot
    

def find_mass(N, N_b, passes, shell_passes, L):
    L_h = 0.1 # length of headers
    m_tubes = 0.2*L*N
    m_nozzles = 0.025*4
    m_shell = 0.65*(L+L_h)
    V_endcaps = 2*np.pi*((0.07/2)**2)*0.002
    rho_plastic = 1150
    m_endcaps = V_endcaps*rho_plastic
    V_tubesheets = V_endcaps*9/2 - N*(0.009*np.pi*(0.008/2)**2)
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
    # print(m,N,N_b,shell_passes, passes, L)
    return m


def test_all(L_list, N_list, N_b_list, passes_list = [1], shell_passes_list = [1], arrange = "tri",Tmethod = "ntu", M = 1.2, year = 2025):
    design_dict = {}
    for N in N_list:
        for N_b in N_b_list:
            for passes in passes_list:
                for shell_passes in shell_passes_list:
                    for L in L_list:
                        d_parameters = (N, N_b, passes, shell_passes, arrange, L)
                        m = find_mass(N, N_b, passes, shell_passes, L)
                        if m < M:
                            if d_parameters[0] < d_parameters[2]:
                                outcome = "more passes than tubes"
                            else:
                                outcome = solution(year, Tmethod, N, N_b, passes, shell_passes, arrange, L)
                            design_dict.update({d_parameters: outcome})
                        elif m >= M:
                            design_dict.update({d_parameters: "too heavy"})
                        else:
                            raise ValueError("invalid mass")
    return design_dict


def find_best(design_dict):
    design_list = list(design_dict.items())
    best_design = ((0,0,0,0,0,0), (0,0))
    for design in design_list:
        # print(design)
        # search for highest heat transfer
        outcome = design[1][1]
        if outcome == "o":
            pass
        elif outcome > float(best_design[1][1]):
            best_design = design
    return best_design

# L = np.linspace(0.15, 0.25, num = 10).tolist()

L = [0.15,0.2,0.25]
N = np.linspace(1, 30, num = 30).tolist()
N_b = np.linspace(0, 50, num =51).tolist()
# passes = np.linspace(1, 6, num =6).tolist()
passes = [ 1,2]
shell_passes = [1,2]

# design_dict = test_all(L, N, N_b, passes, shell_passes, year=2025, Tmethod= "lmtd")
# # print(design_dict)
# best = find_best(design_dict)
# print("------------LMTD method-------------", "\nN:",best[0][0], "\nN_b:",best[0][1], "\npasses:", best[0][2], "\nshell_passes:", best[0][3],  "\narrange:", best[0][4], "\nL:", best[0][5], "\neff:", float(best[1][0]), "\nQdot:",float(best[1][1]) )
# design_dict = test_all(L, N, N_b, passes, shell_passes, year=2025, Tmethod= "ntu")
# best = find_best(design_dict)
# print("------------E-NTU method------------", "\nN:",best[0][0], "\nN_b:",best[0][1], "\npasses:", best[0][2], "\nshell_passes:", best[0][3],  "\narrange:", best[0][4], "\nL:", best[0][5], "\neff:", float(best[1][0]), "\nQdot:",float(best[1][1]) )

# print(solution(2025, "lmtd", 1, 1, 1, 1,'tri', 0.25 ))
# print(solution(2025, "ntu", 16, 6, 2, 2,'tri', 0.2 ))
# print(solution(2025, "ntu", 12, 8, 2, 2,'tri', 0.25 ))
N_b_graph = []
for i in range(len(N_b)):
    sol = solution(2025, "ntu", 12, N_b[i],2,2,'tri',0.25)[1]
    print(sol)
    N_b_graph += [sol]

plt.plot(N_b,N_b_graph)
plt.show()
