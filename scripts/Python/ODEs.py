def ODEs(t, y, Po_ip3r, IcaPQ, volume_ratios, AP_type, stimulus_strength, ISI, N_vgcc,
              coupling_condition, cell_condition, IP3_params, release_params
              ):
    # Solution to ODE's  at previous time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Cac = y[0]  # cytosolic calcium conc
    Ca_ip3 = y[1]  # Nanodomain Calcium
    Ct = y[2]  # ER calcium
    V = y[3]  # Membrane potential
    n2 = y[4]  # Current Gating variable activating
    h2 = y[5]  # Current Gating variable inactivating
    PLC = y[6]  # Extracellular Potassium
    G = y[7]  # Intracellular Sodium
    IP3 = y[8]  # IP3 concentration
    Glut = y[9]  # Glutamate Concentration
    Ca_vgcc = y[10]  # VGCC nanodomain

    ''' Variables for Allosteric Model '''
    # Vesicle number in: Reserve pool - R_allo, Docked pool - U_allo,
    # Vesicle number in: Slow Release Pool (SRP) Refractory pool - RFv_allo, Flow Release Pool (FRP) Refractory pool - RFv_allo
    # Vesicle number in: V(0-1) - SRP, W(0-1) - FRP

    R_allo, U_allo, RFv_allo, RFw_allo = y[11], y[12], y[13], y[14]

    V0, V1, V2, V3, V4, V5 = y[15], y[16], y[17], y[18], y[19], y[20]

    W0, W1, W2, W3, W4, W5 = y[21], y[22], y[23], y[24], y[25], y[26]

    ''' Variables forDual Sensor Model '''
    # Vesicle number in: Reserve pool - R_dual, Docked pool - U_dual,
    # Vesicle number in: Slow Release Pool (SRP) Refractory pool - RFv_dual, Flow Release Pool (FRP) Refractory pool - RFv_dual
    # Vesicle number in: V(0,0 - 5,2) - SRP, W(0,0 - 5,2) - FRP,

    R_dual, U_dual, RFv_dual, RFw_dual = y[27], y[28], y[29], y[30]

    V00, V01, V02, V10, V11, V12, V20, V21, V22 = y[31], y[32], y[33], y[34], y[35], y[36], y[37], y[38], y[39]
    V30, V31, V32, V40, V41, V42, V50, V51, V52 = y[40], y[41], y[42], y[43], y[44], y[45], y[46], y[47], y[48]

    W00, W01, W02, W10, W11, W12, W20, W21, W22 = y[49], y[50], y[51], y[52], y[53], y[54], y[55], y[56], y[57]
    W30, W31, W32, W40, W41, W42, W50, W51, W52 = y[58], y[59], y[60], y[61], y[62], y[63], y[64], y[65], y[66]

    # Parameter Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gamma_1 = volume_ratios["gamma_1"]  # cytoplasmic-to-ER-microdomain Volume ratio
    gamma_2 = volume_ratios["gamma_2"]  # cytoplasmic-to-ER Volume ratio
    gamma_3 = volume_ratios["gamma_3"]  # cytoplasmic-to-VGCC-nanodomain ratio
    Vol_er = 3.9 * 0.1 * 0.1e-18  # (um) ^ 3 ER volume '''not used'''

    ''' Algorithm for applied stimulation '''

    #  Applied/Injected Current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #  Note that the amplitude of applied current is what sets the frequency of the Action Potential (AP) in the train

    if AP_type.lower() == "single":
        # ''' Apply stimulation only after 3 ms. This is chosen so that the system relaxes before stimulation is applied '''
        if t > 3:
            Iapp = 0.0  # uA/cm^2
        else:
            Iapp = stimulus_strength  # uA/cm^2
    elif AP_type.lower() == "train":
        ''' Apply stimulation only for roughly about 400 ms '''
        if t < 3 or t > 435:
            Iapp = 0.0  # uA/cm^2
        else:
            Iapp = stimulus_strength  # uA/cm^2
    elif AP_type.lower() == "paired pulse":
        ''' Employs specified Inter-Spike Interval (ISI) for paired pulse protocol. '''
        ''' Default ISI = 40 ms'''
        if t < 3 or t > 3 + ISI:
            if t > 3 + 3 + ISI:
                Iapp = 0.0  # uA/cm^2
            else:
                Iapp = stimulus_strength
        else:
            Iapp = 0  # uA/cm^2

    # # Nernst potential for Na, K, Cl and Ca %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ENa = 50  # mV
    EK = -100  # mV
    ECl = -70  # mV

    # Currents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ''' Currens From Voltage Gated Channels '''

    gNa = 120  # mS/cm^2  Sodium conductane
    gNaLeak = 0.0175  # mS/cm^2  Sodium Leak conductane
    gK = 36  # mS/cm^2  Potassium Leak conductane
    gKLeak = 0.05  # mS/cm^2  Potassium Leak conductane
    gClLeak = 0.05  # mS/cm^2  Chlorine Leak conductane
    phi = 5.0
    gAHP = 0.01

    # Gating Variables For Current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alphan = 0.01 * (V + 34) / (1 - np.exp(-0.1 * (V + 34)))  # /ms
    betan = 0.125 * np.exp(-(V + 44) / 80)  # /ms
    alphah = 0.07 * np.exp(-(V + 44) / 20)  # /ms
    betah = 1.0 / (1 + np.exp(-0.1 * (V + 14)))  # /ms
    alpham = 0.1 * (V + 30) / (1 - np.exp(-0.1 * (V + 30)))  # /ms
    betam = 4 * np.exp(-(V + 55) / 18)  # /ms
    m1 = alpham * V / (alpham * V + betam * V)

    INa = -gNa * (m1 ** 3) * h2 * (V - ENa) - gNaLeak * (V - ENa)  # uA/cm^2
    IK = -(gK * n2 ** 4 + (gAHP * Cac) / (1 + Cac)) * (V - EK) - gKLeak * (V - EK)  # uA/cm^2
    ICl = -gClLeak * (V - ECl)  # uA/cm^2

    ''' Calculating Endoplasmic Reticulum Calcium '''

    Cer = gamma_2 * (Ct - Cac - Ca_ip3 / gamma_1 - Ca_vgcc / gamma_3)

    ''' Fluxes '''

    ''' PMCA FLUX ####################################################### '''

    Vpmca = 3.195  # uM/ms   Maximum capacity of plasma pump
    Kpmca = 0.5  # uM      Half-maximal activating Cac of plasma pump
    nP = 2.0  # Hill coefficient of plasma pump
    Jpmca = Vpmca * (Cac ** nP) / (Kpmca ** nP + Cac ** nP)  # uM/ms  Flux through plasma pump

    ''' SERCA FLUX ###################################################### '''

    Vs = 10  # uM/ms   Maximum capacity of SERCA                        #0.033;0.9; #2.2;
    Ks = 0.26  # uM    SERCA half-maximal activating Cac                #0.45;0.5; #0.1;
    ns = 1.75  # Hill coefficient of SERCA

    Jserca = Vs * (Cac ** ns) / (Ks ** ns + Cac ** ns)  # uM/ms  Uptake of Ca into the ER by SERCA pumps

    ''' VGCC FLUX ####################################################### '''
    Farad = 96485.33  # C/mole  Faraday's constant, unit: coul/mole
    z = 2  # valence of Ca ion
    Vol_tmnal = 1.22e-16  # Litre    Presynaptic terminal Volume (Assuming a spherical Bouton shape)
    cluster_radius = 25e-03  # um
    cluster_area = (np.pi()) * (cluster_radius) ^ 2  # um^2
    active_zone_area = 0.04  # um^2
    active_zone_number = 1.3  # Number of active zones
    Area = 3.8489e-09  # cm^2   Bouton Membrane area
    Vol_tmnal = 1.7962e-16  # Litre  Bouton Volume (Assuming a spherical Bouton shape)
    Kdiff_vgcc = 0.071  # /ms  Rate of diffusion from VGCC nanodomain

    channel_density = N_vgcc / (active_zone_number * active_zone_area)
    Ica = channel_density * cluster_area * IcaPQ

    #  Macroscopic current density through an ensemble of open channels where channel_density is crucial:
    Jvgcc = (-Ica / (
            z * Farad * Vol_tmnal)) * 1e-03  # uM/ms  Calcium influx from the extracellular space to cytosol(with Mitoc 2.03)
    Jdiff_vgcc = Kdiff_vgcc * (Ca_vgcc - Cac)  # uM/ms  Diffusion from VGCC cluster nanodomain to the cytoplasm

    ''' RECEPTOR FLUXES --- IP3 and RYR ################################# '''

    KIPR = 5  # /ms     IP3R flux coefficient
    kdiff = 10  # /ms     Ca diffusional flux coefficient

    Jipr =  KIPR * Po_ip3r * (Cer - Ca_ip3)  # uM/ms  Flux through the IP3R
    Jdiff = kdiff * (Ca_ip3 - Cac)  # uM/ms  Diffusion from ER cluster nanodomain to the cytoplasm

    """ Leak Fluxes -- ER and Plasma Membrane ############################# """

    kleak = 0.0022  # /ms     ER leak flux coefficient                         #0.0032*1e-03; #0.0032 #0.1221;
    Jleakin = 0.03115  # uM/ms   Plasma membrane leak influx
    Vleakin = 0.2  # /ms     IP3 In Leak flux coefficient

    Jleak = kleak * (Cer - Cac)  # uM/ms  Baground leak from the ER into the cytoplasm
    Jin = Jleakin + Vleakin * IP3  # uM/ms  Sum of main Ca2+ influxes

    """ Flux from VGCC Nano Domain to IP3R Nanodomain ###################### """

    if coupling_condition == "Higher_Coupling_AD":
        if cell_condition.lower() == "wt":
            Kc = 20
            Vc = 118
            K_hat = 5
        elif cell_condition == "ad":
            Kc = 10
            Vc = 118
            K_hat = 15
    elif coupling_condition == "Higher_Coupling_WT":
        if cell_condition == "wt":
            Kc = 10
            Vc = 118
            K_hat = 15
        elif cell_condition == "ad":
            Kc = 20
            Vc = 118
            K_hat = 5
    elif coupling_condition == "Same_Coupling":
            Kc = 10
            Vc = 118
            K_hat = 15

    JCa_vgcc_Ca_ip3 = Vc * (Ca_vgcc ** 2 - K_hat * Ca_ip3 ** 2) / (Kc ** 2 + Ca_vgcc ** 2)

    ''' Differential Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '''

    """ Membrane Potential ########################################### """

    C = 1  # uF/cm^2  membrane capacitance
    dvdt = (Iapp + INa + IK + ICl) / C + (Jin * z * Farad * Vol_tmnal) / (C * Area)  # Jpm and Jin times factor z*F

    """ Gating Variables ############################################# """

    dn2dt = phi * (alphan * (1 - n2) - betan * n2)
    dh2dt = phi * (alphah * (1 - h2) - betah * h2)

    """ Calcium Dynamics ############################################## """

    dCacdt = Jdiff + Jdiff_vgcc + Jleak + Jin - Jpmca - Jserca
    dCa_ip3dt = gamma_1 * (Jipr - Jdiff) + JCa_vgcc_Ca_ip3
    dCa_vgccdt = gamma_3 * (Jvgcc - Jdiff_vgcc - JCa_vgcc_Ca_ip3 / gamma_1)
    dCtdt = (Jin + Jvgcc - Jpmca)

    """ IP3 Production ################################################ """

    # Parameters from Joe Latulippe. et. al 2021
    # Contains Amyloid beta effect.

    G_tot = 1  # Scaled total number of G-protein
    PLC_tot = 1  # Scaled total number of PLC
    abeta = 0  # A? concentration
    rbeta = 0.001  # A? decay rate
    Kq = 0.0086  # ug/mL PLC dissociation constant

    [Vo, Vq, Kip3k, Kplc, k_3k, k_5p,
     Kf_plc, Kb_plc, Kf_gp, Kb_gp, delta_G, Vr, Kr] = IP3_params[
        cell_condition.lower()].values()

    t_ip3 = 1 / (k_3k + k_5p)  # ms    IP3 represents the characteristic time of IP3 turnover
    Nu = k_3k / (k_3k + k_5p)  # affects negative feedback from Ca2+ to IP3 degredation
    q = np.heaviside(t - 2, 0.5) * abeta * np.exp(-rbeta * (t - 2) * np.heaviside(t - 2, 0.5))
    rho_G = Vr * q / (Kr + q)  # A? dependent G-protein activation
    Vplc = Vo + (Vq * q ** 2 / (Kq ** 2 + q ** 2))  # A? dependent maximal rate of PLC mediated IP3 production
    Jplc = Vplc * PLC * (Cac ** 2 / (Kplc ** 2 + Cac ** 2))
    Jdeg = (Nu * (Cac ** 2 / (Kip3k ** 2 + Cac ** 2)) + (1 - Nu)) * IP3
    dPLCdt = Kf_plc * G * (PLC_tot - PLC) - Kb_plc * PLC
    dGdt = Kf_gp * (rho_G + delta_G) * (G_tot - G) - Kb_gp * G
    dIP3dt = 0 * (Jplc - Jdeg) / t_ip3
    dGlutdt = 0

    """ Transmitter Release Rate Equations """

    """ %% Allosteric Model Paramters """

    [Kon_allo, Koff_allo, I_allo, F_allo, b_allo, Krf_allo,
     Kmob_allo, Kdemob_allo, Kprime_allo, Kupr_allo, Kattach_allo, Kdetach_allo] = release_params[
        'allo'].values()  # /ms

    # Allosteric Model

    # Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Vtotal_allo = V0 + V1 + V2 + V3 + V4 + V5

    Wtotal_allo = W0 + W1 + W2 + W3 + W4 + W5

    dR_allodt = -Kmob_allo * Cac * R_allo + U_allo * Kdemob_allo
    dU_allodt = -Kdemob_allo * U_allo + Kmob_allo * Cac * R_allo - Kprime_allo * U_allo * Cac * (
            1 - RFv_allo) + Kupr_allo * V0

    dRFv_allodt = (1 - RFv_allo) * (I_allo * (V0 + V1 * F_allo + V2 * F_allo ** 2 + V3 *
                                              F_allo ** 3 + V4 * F_allo ** 4 + V5 *
                                              F_allo ** 5)) / Vtotal_allo - Krf_allo * RFv_allo

    dRFw_allodt = (1 - RFw_allo) * (
            I_allo * (
            W0 + W1 * F_allo + W2 * F_allo ** 2 + W3 *
            F_allo ** 3 + W4 * F_allo ** 4 + W5 *
            F_allo ** 5)) / Wtotal_allo - Krf_allo * RFw_allo

    dV0dt = Kprime_allo * U_allo * Cac * (1 - RFv_allo) - Kupr_allo * V0 - \
            Kattach_allo * V0 * Ca_vgcc * (1 - RFw_allo) + Kdetach_allo * W0 + \
            Koff_allo * V1 \
            - I_allo * V0 - 5 * Kon_allo * Cac * V0
    dV1dt = 2 * Koff_allo * b_allo * V2 - Kon_allo * V1 + 5 * Kon_allo * Cac * V0 - 4 * Kon_allo * Cac * V1 - I_allo * F_allo * V1
    dV2dt = 3 * Koff_allo * V3 * b_allo ** 2 - 2 * \
            Koff_allo * V2 * b_allo + 4 * Kon_allo * Cac * V1 - 3 * Kon_allo * Cac * V2 - I_allo * V2 * F_allo ** 2
    dV3dt = 4 * Koff_allo * V4 * b_allo ** 3 - 3 * \
            Koff_allo * V3 * b_allo ** 2 + 3 * Kon_allo * Cac * V2 - 2 * Kon_allo * Cac * V3 - I_allo * V3 * F_allo ** 3
    dV4dt = 5 * Koff_allo * V5 * b_allo ** 4 - 4 * \
            Koff_allo * V4 * b_allo ** 3 + 2 * Kon_allo * Cac * V3 - Kon_allo * Cac * V4 - I_allo * V4 * \
            F_allo ** 4
    dV5dt = Kon_allo * Cac * V4 - 5 * Koff_allo * V5 * \
            b_allo ** 4 - I_allo * V5 * F_allo ** 5

    dW0dt = Kattach_allo * W0 * Ca_vgcc * (
            1 - RFw_allo) - Kdetach_allo * W0 + Koff_allo * W1 - I_allo * W0 - 5 * Kon_allo * Ca_vgcc * W0
    dW1dt = 2 * Koff_allo * b_allo * W2 - Kon_allo * W1 + 5 * Kon_allo * Ca_vgcc * W0 - 4 * Kon_allo * Ca_vgcc * W1 - I_allo * F_allo * W1
    dW2dt = 3 * Koff_allo * W3 * b_allo ** 2 - 2 * \
            Koff_allo * W2 * b_allo + 4 * Kon_allo * Ca_vgcc * W1 - 3 * Kon_allo * Ca_vgcc * W2 - I_allo * W2 * F_allo ** 2
    dW3dt = 4 * Koff_allo * W4 * b_allo ** 3 - 3 * \
            Koff_allo * W3 * b_allo ** 2 + 3 * Kon_allo * Ca_vgcc * W2 - 2 * Kon_allo * Ca_vgcc * W3 - I_allo * W3 * F_allo ** 3
    dW4dt = 5 * Koff_allo * W5 * b_allo ** 4 - 4 * \
            Koff_allo * W4 * b_allo ** 3 + 2 * Kon_allo * Ca_vgcc * W3 - Kon_allo * Ca_vgcc * W4 - I_allo * W4 * \
            F_allo ** 4
    dW5dt = Kon_allo * Ca_vgcc * W4 - 5 * Koff_allo * W5 * \
            b_allo ** 4 - I_allo * W5 * F_allo ** 5

    """ Dual Censor Model  ################################### """

    # Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [alpha_dual, beta_dual, Chi_dual, delta_dual, a_dual, bo_dual, gam_2_dual, gam_3_dual, gam_1_dual, Krf_dual,
     Kmob_dual, Kdemob_dual, Kprime_dual, Kupr_dual, Kattach_dual, Kdetach_dual] = release_params['dual'].values()

    Vtotal_dual = V00 + V01 + V02 + V10 + V11 + V12 + V20 + V21 + V22 + V30 + V31 + V32 + V40 + V41 + V42 + V50 + V51 + V52

    Wtotal_dual = W00 + W01 + W02 + W10 + W11 + W12 + W20 + W21 + W22 + W30 + W31 + W32 + W40 + W41 + W42 + W50 + W51 + W52

    dR_dualdt = -Kmob_dual * Cac * R_dual + U_dual * Kdemob_dual
    dU_dualdt = -Kprime_dual * U_dual * Cac * (1 - RFv_dual) + Kupr_dual * V00

    dRFv_dualdt = (1 - RFv_dual) * (
            gam_3_dual * (V02 + V12 + V22 + V32 + V42 + V52) + gam_2_dual * (V50 + V51 + V52) + gam_1_dual * V00) / \
                  Vtotal_dual - Krf_dual * RFv_dual

    dRFw_dualdt = (1 - RFw_dual) * (
            gam_3_dual * (W02 + W12 + W22 + W32 + W42 + W52) + gam_2_dual * (W50 + W51 + W52) + gam_1_dual * W00) / \
                  Wtotal_dual - Krf_dual * RFw_dual

    dV00 = Kprime_dual * U_dual * Cac * (1 - RFv_dual) - Kupr_dual * V00 - \
           Kattach_dual * V00 * Ca_vgcc * (
                   1 - RFw_dual) + Kdetach_dual * W00 + beta_dual * \
           V10 - 5 * alpha_dual * V00 * Cac + delta_dual * V01 - 2 * \
           Chi_dual * V00 * Cac - gam_1_dual * V00
    dV01 = 2 * Chi_dual * V00 * Cac - delta_dual * V01 + 2 * \
           bo_dual * delta_dual * V02 - Chi_dual * Cac * V01 + beta_dual * V11 - 5 * alpha_dual * Cac * V01
    dV02 = Chi_dual * Cac * V01 - 2 * bo_dual * delta_dual * V02 - gam_3_dual * V02 - 5 * alpha_dual * Cac * V02 + \
           beta_dual * V12
    dV10 = 5 * alpha_dual * Cac * V00 - beta_dual * V10 - 4 * \
           alpha_dual * Cac * V10 + 2 * bo_dual * \
           beta_dual * V20 - 2 * Chi_dual * Cac * V10 + \
           delta_dual * V11
    dV11 = 5 * alpha_dual * Cac * V01 - beta_dual * V11 - 4 * \
           alpha_dual * Cac * V11 + 2 * bo_dual * \
           beta_dual * V21 + 2 * Chi_dual * Cac * V10 - \
           delta_dual * \
           V11 - Chi_dual * Cac * V11 + 2 * bo_dual * \
           delta_dual * V12
    dV12 = 5 * alpha_dual * Cac * V02 - beta_dual * V12 - 4 * \
           alpha_dual * Cac * V12 + 2 * bo_dual * \
           beta_dual * V22 + Chi_dual * Cac * V11 - 2 * \
           bo_dual * \
           delta_dual * V12 - gam_3_dual * V12
    dV20 = 4 * alpha_dual * Cac * V10 - 2 * bo_dual * \
           beta_dual * V20 - 3 * alpha_dual * Cac * V20 + 3 * (
                   bo_dual ** 2) * beta_dual * V30 - 2 * \
           Chi_dual * \
           Cac * V20 + delta_dual * V21
    dV21 = 4 * alpha_dual * Cac * V11 - 2 * bo_dual * \
           beta_dual * V21 - 3 * alpha_dual * Cac * V21 + 3 * (
                   bo_dual ** 2) * beta_dual * V31 + 2 * \
           Chi_dual * \
           Cac * V20 - delta_dual * V21 - Chi_dual * Cac * V21 + 2 * \
           bo_dual * delta_dual * V22
    dV22 = 4 * alpha_dual * Cac * V12 - 2 * bo_dual * \
           beta_dual * V22 - 3 * alpha_dual * Cac * V22 + 3 * (
                   bo_dual ** 2) * beta_dual * V32 + \
           Chi_dual * Cac * \
           V21 - 2 * bo_dual * delta_dual * V22 - \
           gam_3_dual * V22
    dV30 = 3 * alpha_dual * Cac * V20 - 3 * (bo_dual ** 2) * \
           beta_dual * V30 - 2 * alpha_dual * Cac * V30 + 4 * (
                   bo_dual ** 3) * beta_dual * V40 - 2 * \
           Chi_dual * Cac * V30 + delta_dual * V31
    dV31 = 3 * alpha_dual * Cac * V21 - 3 * (bo_dual ** 2) * \
           beta_dual * V31 - 2 * alpha_dual * Cac * V31 + 4 * (
                   bo_dual ** 3) * beta_dual * V41 + 2 * \
           Chi_dual * Cac * V30 - delta_dual * V31 - \
           Chi_dual * Cac * V31 + 2 * bo_dual * delta_dual * V32
    dV32 = 3 * alpha_dual * Cac * V22 - 3 * (bo_dual ** 2) * \
           beta_dual * V32 - 2 * alpha_dual * Cac * V32 + 4 * (
                   bo_dual ** 2) * beta_dual * V42 + \
           Chi_dual * \
           Cac * V31 - 2 * bo_dual * delta_dual * V32 - \
           gam_3_dual * V32
    dV40 = 2 * alpha_dual * Cac * V30 - 4 * (bo_dual ** 3) * \
           beta_dual * V40 - alpha_dual * Cac * V40 + 5 * (
                   bo_dual ** 4) * beta_dual * V50 - 2 * \
           Chi_dual * \
           Cac * V40 + delta_dual * V41
    dV41 = 2 * alpha_dual * Cac * V31 - 4 * (bo_dual ** 3) * \
           beta_dual * V41 - alpha_dual * Cac * V41 + 5 * (
                   bo_dual ** 4) * beta_dual * V51 + 2 * \
           Chi_dual * \
           Cac * V40 - delta_dual * V41 - Chi_dual * Cac * V41 + 2 * \
           bo_dual * delta_dual * V42
    dV42 = 2 * alpha_dual * Cac * V32 - 4 * (bo_dual ** 3) * \
           beta_dual * V42 - alpha_dual * Cac * V42 + 5 * (
                   bo_dual ** 4) * beta_dual * V52 + \
           Chi_dual * Cac * \
           V41 - 2 * bo_dual * delta_dual * V42 - \
           gam_3_dual * V42
    dV50 = alpha_dual * Cac * V40 - 5 * (bo_dual ** 4) * \
           beta_dual * V50 - 2 * Chi_dual * Cac * V50 + \
           delta_dual * V51 - gam_2_dual * V50
    dV51 = alpha_dual * Cac * V41 - 5 * (
            bo_dual ** 4) * beta_dual * V51 + 2 * \
           Chi_dual * Cac * V50 - delta_dual * V51 - \
           Chi_dual * Cac * V51 + 2 * bo_dual * \
           delta_dual * V52 - gam_2_dual * V51
    dV52 = alpha_dual * Cac * V42 - 5 * (
            bo_dual ** 4) * beta_dual * V52 + Chi_dual * Cac * V51 - 2 * bo_dual * delta_dual * V52 - \
           gam_3_dual * V52 - gam_2_dual * V52

    dW00 = Kattach_dual * V00 * Ca_vgcc * (
            1 - RFw_dual) - Kdetach_dual * W00 + beta_dual * W10 - 5 * \
           alpha_dual * W00 * Ca_vgcc + delta_dual * W01 - 2 * \
           Chi_dual * \
           W00 * Ca_vgcc - gam_1_dual * W00
    dW01 = 2 * Chi_dual * W00 * Ca_vgcc - delta_dual * W01 + 2 * \
           bo_dual * delta_dual * W02 - Chi_dual * Ca_vgcc * W01 + beta_dual * W11 - 5 * alpha_dual * Ca_vgcc * W01
    dW02 = Chi_dual * Ca_vgcc * W01 - 2 * bo_dual * delta_dual * W02 - gam_3_dual * W02 - 5 * alpha_dual * Ca_vgcc * W02 + \
           beta_dual * W12
    dW10 = 5 * alpha_dual * Ca_vgcc * W00 - beta_dual * W10 - 4 * \
           alpha_dual * Ca_vgcc * W10 + 2 * bo_dual * \
           beta_dual * W20 - 2 * Chi_dual * Ca_vgcc * W10 + \
           delta_dual * W11
    dW11 = 5 * alpha_dual * Ca_vgcc * W01 - beta_dual * W11 - 4 * \
           alpha_dual * Ca_vgcc * W11 + 2 * bo_dual * \
           beta_dual * W21 + 2 * Chi_dual * Ca_vgcc * W10 - \
           delta_dual * \
           W11 - Chi_dual * Ca_vgcc * W11 + 2 * bo_dual * \
           delta_dual * W12
    dW12 = 5 * alpha_dual * Ca_vgcc * W02 - beta_dual * W12 - 4 * \
           alpha_dual * Ca_vgcc * W12 + 2 * bo_dual * \
           beta_dual * W22 + Chi_dual * Ca_vgcc * W11 - 2 * \
           bo_dual * \
           delta_dual * W12 - gam_3_dual * W12
    dW20 = 4 * alpha_dual * Ca_vgcc * W10 - 2 * bo_dual * \
           beta_dual * W20 - 3 * alpha_dual * Ca_vgcc * W20 + 3 * (
                   bo_dual ** 2) * beta_dual * W30 - 2 * \
           Chi_dual * \
           Ca_vgcc * W20 + delta_dual * W21
    dW21 = 4 * alpha_dual * Ca_vgcc * W11 - 2 * bo_dual * \
           beta_dual * W21 - 3 * alpha_dual * Ca_vgcc * W21 + 3 * (
                   bo_dual ** 2) * beta_dual * W31 + 2 * \
           Chi_dual * \
           Ca_vgcc * W20 - delta_dual * W21 - Chi_dual * Ca_vgcc * W21 + 2 * \
           bo_dual * delta_dual * W22
    dW22 = 4 * alpha_dual * Ca_vgcc * W12 - 2 * bo_dual * \
           beta_dual * W22 - 3 * alpha_dual * Ca_vgcc * W22 + 3 * (
                   bo_dual ** 2) * beta_dual * W32 + \
           Chi_dual * Ca_vgcc * \
           W21 - 2 * bo_dual * delta_dual * W22 - \
           gam_3_dual * W22
    dW30 = 3 * alpha_dual * Ca_vgcc * W20 - 3 * (bo_dual ** 2) * \
           beta_dual * W30 - 2 * alpha_dual * Ca_vgcc * W30 + 4 * (
                   bo_dual ** 3) * beta_dual * W40 - 2 * \
           Chi_dual * Ca_vgcc * W30 + delta_dual * W31
    dW31 = 3 * alpha_dual * Ca_vgcc * W21 - 3 * (bo_dual ** 2) * \
           beta_dual * W31 - 2 * alpha_dual * Ca_vgcc * W31 + 4 * (
                   bo_dual ** 3) * beta_dual * W41 + 2 * \
           Chi_dual * Ca_vgcc * W30 - delta_dual * W31 - \
           Chi_dual * Ca_vgcc * W31 + 2 * bo_dual * delta_dual * W32
    dW32 = 3 * alpha_dual * Ca_vgcc * W22 - 3 * (bo_dual ** 2) * \
           beta_dual * W32 - 2 * alpha_dual * Ca_vgcc * W32 + 4 * (
                   bo_dual ** 2) * beta_dual * W42 + \
           Chi_dual * \
           Ca_vgcc * W31 - 2 * bo_dual * delta_dual * W32 - \
           gam_3_dual * W32
    dW40 = 2 * alpha_dual * Ca_vgcc * W30 - 4 * (bo_dual ** 3) * \
           beta_dual * W40 - alpha_dual * Ca_vgcc * W40 + 5 * (
                   bo_dual ** 4) * beta_dual * W50 - 2 * \
           Chi_dual * \
           Ca_vgcc * W40 + delta_dual * W41
    dW41 = 2 * alpha_dual * Ca_vgcc * W31 - 4 * (
            bo_dual ** 3) * beta_dual * W41 - alpha_dual * Ca_vgcc * W41 + 5 * (
                   bo_dual ** 4) * beta_dual * W51 + 2 * \
           Chi_dual * \
           Ca_vgcc * W40 - delta_dual * W41 - Chi_dual * Ca_vgcc * W41 + 2 * \
           bo_dual * delta_dual * W42
    dW42 = 2 * alpha_dual * Ca_vgcc * W32 - 4 * (bo_dual ** 3) * \
           beta_dual * W42 - alpha_dual * Ca_vgcc * W42 + 5 * (
                   bo_dual ** 4) * beta_dual * W52 + \
           Chi_dual * Ca_vgcc \
           * W41 - 2 * bo_dual * delta_dual * W42 - \
           gam_3_dual * W42
    dW50 = alpha_dual * Ca_vgcc * W40 - 5 * (bo_dual ** 4) * \
           beta_dual * W50 - 2 * Chi_dual * Ca_vgcc * W50 + \
           delta_dual * W51 - gam_2_dual * W50
    dW51 = alpha_dual * Ca_vgcc * W41 - 5 * (
            bo_dual ** 4) * beta_dual * W51 + 2 * \
           Chi_dual * Ca_vgcc * W50 - delta_dual * W51 - \
           Chi_dual * Ca_vgcc * W51 + 2 * bo_dual * \
           delta_dual * W52 - gam_2_dual * W51
    dW52 = alpha_dual * Ca_vgcc * W42 - 5 * (
            bo_dual ** 4) * beta_dual * W52 + Chi_dual * Ca_vgcc * W51 - 2 * bo_dual * delta_dual * W52 - \
           gam_3_dual * W52 - gam_2_dual * W52

    # # Return Solutions

    # Deterministic Solutions
    dydt = [dCacdt, dCa_ip3dt, dCtdt, dvdt, dn2dt, dh2dt,
            dPLCdt, dGdt, dIP3dt, dGlutdt,
            dCa_vgccdt,
            dR_allodt, dU_allodt, dRFv_allodt, dRFw_allodt,
            dV0dt, dV1dt, dV2dt, dV3dt, dV4dt, dV5dt, dW0dt, dW1dt, dW2dt, dW3dt, dW4dt, dW5dt,
            dR_dualdt, dU_dualdt, dRFv_dualdt, dRFw_dualdt,
            dV00, dV01, dV02, dV10, dV11, dV12, dV20, dV21, dV22, dV30, dV31, dV32, dV40, dV41, dV42, dV50, dV51,
            dV52,
            dW00, dW01, dW02, dW10, dW11, dW12, dW20, dW21, dW22, dW30, dW31, dW32, dW40, dW41, dW42, dW50, dW51,
            dW52]

    return dydt