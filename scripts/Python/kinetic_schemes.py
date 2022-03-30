
def kinetic_schemes(t, y, dt, state_ip3r, state_PQ, N_vgcc, N_ipr, IP3Receptor_params, cell_condition):
    Ca_ip3 = y[1]  # Nanodomain Calcium
    Ca_vgcc = y[10]
    V = y[3]  # Membrane potential
    IP3 = y[8]  # IP3 Concentration
    Cout = 2000  # Extracellular Ca concentration, unit: uM
    Farad = 96485.33  # Coul/mole         Faraday's constant, unit: coul/mole
    T = 300  # K                 Temperature, unit:Kelvin
    Ro = 8.31  # J/(mole*K)        Ideal gas constant, unit: J/(mole*K)
    z = 2  # valence of Ca ion

    # # Nernst potential for Ca %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Eca = (Ro * T / (z * Farad)) * 1000 * np.log(Cout / Ca_vgcc)  # mV

    # %%%%%%%%%%%%%%  IP3 Kinetic Scheme %%%%%%%%%%%%%%%%%%%%%

    [a1, nO, Kod, a2, nA, Kad, a3, nI,
     Kid, j01, j12, j22, j23, j45, J01_tilda, J45_tilda] = IP3Receptor_params[cell_condition.lower()].values()

    # # Transition Rates and Open Probability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    KO = (a1 * IP3 ** nO) / (IP3 ** nO + Kod ** nO)
    KI = (a3 * IP3 ** nI) / (IP3 ** nI + Kid ** nI)
    KA = (a2 * IP3 ** nA) / (IP3 ** nA + Kad ** nA)

    kRA = 1 / ((1 / (j01 * Ca_ip3)) + (1 / (j12 * (Ca_ip3 ** 2))))
    kAR = 1 / (KA * (Ca_ip3 ** 2) * ((1 / (j01 * Ca_ip3)) + (1 / (j12 * (Ca_ip3 ** 2)))))

    kAO = 1 / (KA * (Ca_ip3 ** 2) * (1 / (j22 * (Ca_ip3 ** 2))))
    kOA = 1 / (KO * (Ca_ip3 ** 2) * (1 / (j22 * (Ca_ip3 ** 2))))

    kOI = 1 / ((KO * Ca_ip3 ** 2) * (1 / (j23 * Ca_ip3 ** 3) + 1 / (j45 * Ca_ip3 ** 5)))
    kIO = 1 / ((KI * Ca_ip3 ** 5) * (1 / (j23 * Ca_ip3 ** 3) + 1 / (j45 * Ca_ip3 ** 5)))

    kRI = 1 / (1 / (J01_tilda * Ca_ip3) + 1 / (J45_tilda * Ca_ip3 ** 5))
    kIR = 1 / ((KI * Ca_ip3 ** 5) * (1 / (J01_tilda * Ca_ip3) + 1 / (J45_tilda * Ca_ip3 ** 5)))

    # # For Stochastic IP3R Channel Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pRA = kRA * dt
    pAR = kAR * dt
    pOA = kOA * dt
    pAO = kAO * dt
    pRI = kRI * dt
    pIR = kIR * dt
    pIO = kIO * dt
    pOI = kOI * dt

    OpenNPR = 0
    for i in range(0, N_ipr):
        z = random()
        if state_ip3r[0, i] == 1:
            if z < pRA:
                state_ip3r[0, i] = 2
            elif pRA <= z < (pRA + pRI):
                state_ip3r[0, i] = 4
        elif state_ip3r[0, i] == 2:
            if z < pAR:
                state_ip3r[0, i] = 1
            elif pAR <= z < (pAR + pAO):
                state_ip3r[0, i] = 3
        elif state_ip3r[0, i] == 3:
            if z < pOA:
                state_ip3r[0, i] = 2
            elif pOA <= z < (pOA + pOI):
                state_ip3r[0, i] = 4
        elif state_ip3r[0, i] == 4:
            if z < pIR:
                state_ip3r[0, i] = 1
            elif pIR <= z < (pIR + pIO):
                state_ip3r[0, i] = 3

        if state_ip3r[0, i] == 3:
            OpenNPR = OpenNPR + 1

    Po_ip3r = (OpenNPR / N_ipr)

    # #  Stochastic VGCC channel simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Suhita Nadkarni P-Q- Type current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gPQ = 2.7e-09  # mS

    alpha_10PQ = 4.04  # /ms
    beta_10PQ = 2.88  # /ms
    k1_PQ = 49.14  # mV
    alpha_20PQ = 6.70  # /ms
    beta_20PQ = 6.30  # /ms
    k2_PQ = 42.08  # mV
    alpha_30PQ = 4.39  # /ms
    beta_30PQ = 8.16  # /ms
    k3_PQ = 55.31  # mV
    alpha_40PQ = 17.33  # /ms
    beta_40PQ = 1.84  # /ms
    k4_PQ = 26.55  # mV

    alpha_1PQ = alpha_10PQ * np.exp(V / k1_PQ)  # /ms
    beta_1PQ = beta_10PQ * np.exp(-V / k1_PQ)  # /ms
    alpha_2PQ = alpha_20PQ * np.exp(V / k2_PQ)  # /ms
    beta_2PQ = beta_20PQ * np.exp(-V / k2_PQ)  # /ms
    alpha_3PQ = alpha_30PQ * np.exp(V / k3_PQ)  # /ms
    beta_3PQ = beta_30PQ * np.exp(-V / k3_PQ)  # /ms
    alpha_4PQ = alpha_40PQ * np.exp(V / k4_PQ)  # /ms
    beta_4PQ = beta_40PQ * np.exp(-V / k4_PQ)  # /ms

    OpenPQ = 0
    for i in range(0, N_vgcc):
        z = random()
        if state_PQ[0, i] == 0:
            if z < alpha_1PQ * dt:
                state_PQ[0, i] = 1
        elif state_PQ[0, i] == 1:
            if z < alpha_2PQ * dt:
                state_PQ[0, i] = 2
            elif alpha_2PQ * dt <= z < (alpha_2PQ + beta_1PQ) * dt:
                state_PQ[0, i] = 0
        elif state_PQ[0, i] == 2:
            if z < alpha_3PQ * dt:
                state_PQ[0, i] = 3
            elif alpha_3PQ * dt <= z < (alpha_3PQ + beta_2PQ) * dt:
                state_PQ[0, i] = 1
        elif state_PQ[0, i] == 3:
            if z < alpha_4PQ * dt:
                state_PQ[0, i] = 4
            elif alpha_4PQ * dt <= z < (alpha_4PQ + beta_3PQ) * dt:
                state_PQ[0, i] = 2
        elif state_PQ[0, i] == 4:
            if z < beta_4PQ * dt:
                state_PQ[0, i] = 3

        if state_PQ[0, i] == 4:
            OpenPQ = OpenPQ + 1

    PoPQ = OpenPQ / N_vgcc
    IcaPQ = gPQ * PoPQ * (V - Eca)  # uA

    return [PoPQ, Po_ip3r, IcaPQ, state_ip3r, state_PQ]