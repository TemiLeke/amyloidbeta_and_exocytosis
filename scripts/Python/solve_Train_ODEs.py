import numpy as np
import os
import scipy.integrate
import pandas as pd
import random
from kinetic_schemes import kinetic_schemes
from ODEs import ODEs
from solve_ODEs import solve_ODEs


#        Description

#       This script simulates calcium dynamics, neurotransmitter release in
#       physiologically plausible hippocampal CA1 terminal mode.
#       All simulations have a time step of 1 us (i,e 0.001 ms). We
#       simulated fo 50 trials and calculated average response for VGGC
#       values in the range= 5--150 (represented as k in simulation).
#       Release probability in response to an AP is calculated by
#       integrating the release rate for slow and fast
#       vesicles divided by the number of RRV initially in both pools.
#       Simulations for single AP and paired pulse protocol were run for
#       100ms, while AP Train simulation were done for 450 ms. For paired
#       pulse protocol, the total duration of each pulse is 30ms. ISI can be
#       modified as desired. AP train of 20Hz was used in our simulations.
#       Simulation results (average of 50 trials) are saved to an automatically created data
#       directory for each VGCC number.
#       The simulation runs for all coupling configurations
#       (High Coupling and Low Coupling as describes in paper) and cell conditions (WT and AD).


###### Model Parameters and Configuration

#       All kinetic and reaction rates/parameters are specified in the
#       kinetic_schemes.py and ODEs.py file. Each parameter is annotated
#       accordingly.
#       All state variables are set to closed, and all dynamical variables
#       are set to physiologically plausible values. Current stimulation is
#       only applied 3 ms after initialization of simulation. This ensures
#       that all variables attain their steady-state resting levels.
#       The coupled ODE system is solved using the RK4 algorithm as other
#       inbuilt methods appear unstable with stochastic kinetic schemes.



####### Note :

#       Please note that all files; "kinetic_schemes.py", "ODEs.py", "solve_ODEs.py"
#       "solve_PPR_ODEs.py", "solve_Train_ODEs.py", "solve_singleAP_ODEs.py", and
#       have to be in the same directory for
#       simulation to run without hiccups.


if __name__ == "__main__":

    N_ipr = 10  # number of IP3 receptors
    ISI = 40  # ms  Inter-spike Interval for Paired-Pulse protocol
    nt = 450000  # numbser of time steps
    dT_save = 0.001  # time step for saving data
    time_step = 0.001
    vgcc_channel_number_range = range(5, 155, 5)
    volume_ratios = {"gamma_1": 100, "gamma_2": 10, "gamma_3": 60}
    cell_conditions = ["wt", "ad"]
    coupling_conditions = ["Higher_Coupling_WT", "Higher_Coupling_AD", "Same_Coupling"] #  Coupling configuration
    num_trials = 50


    init = [0.1, 0.1, 56.0, -70, 0.01, 0.01,                      # dCacdt, dCa_iprdt, dCtdt, dvdt, dn2dt, dh2dt
            1, 1, 0.16, 15,                                       # dPLCdt, dGdt, dIP3dt, dGlutdt,
            0.1,                                                  # dCa_vgccdt
            170, 20, 0, 0, 5, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0,    # dR_allodt, dU_allodt, dRFv_allodt, dRFw_allodt, dV0dt, dV1dt, dV2dt, dV3dt, dV4dt, dV5dt, dW0dt, dW1dt, dW2dt, dW3dt, dW4dt, dW5dt
            170, 20, 0, 0,                                        # dR_dualdt, dU_dualdt, dRFv_dualdt, dRFw_dualdt
            5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # dV00, dV01, dV02, dV10, dV11, dV12, dV20, dV21, dV22, dV30, dV31, dV32, dV40, dV41, dV42, dV50, dV51, dV52,
            5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # dW00, dW01, dW02, dW10, dW11, dW12, dW20, dW21, dW22, dW30, dW31, dW32, dW40, dW41, dW42, dW50, dW51, dW52,

    IP3Receptor_params = {"wt": {"a1": 17.050543, "nO": 2.473407, "Kod": 0.909078, "a2": 18.49186, "nA": 0.093452,
                               "Kad": 1.955650e03, "a3": 2.73028e02, "nI": 56.84823, "Kid": 0.089938,
                               "j01": 3.031635e02, "j12": 3.230063e02, "j22": 4.814111, "j23": 5.356155,
                               "j45": 5.625616, "J01_tilda": 3.013284e02, "J45_tilda": 2.648741
                               },
                        "ad": {"a1": 1.108278e02, "nO": 2.473407, "Kod": 0.909078, "a2": 18.49186, "nA": 0.093452,
                                    "Kad": 1.955650e03, "a3": 2.73028e02, "nI": 56.84823, "Kid": 0.089938,
                                    "j01": 3.031635e02, "j12": 3.230063e02, "j22": 5.3978052, "j23": 2.0652269e03,
                                    "j45": 5.4319289, "J01_tilda": 3.013284e02, "J45_tilda": 8.512829e-08
                                    }
                        }
    IP3_params = {"wt": {"Vo": 0.15,  # uM    Intrinsic PLC-mediated IP3 production
                                   "Vq": 7.82,  # uM    Control parameter for influence of A? on IP3
                                   "Kip3k": 0.6,  # uM    Half-activation for 3-kinase
                                   "Kplc": 0.01,  # uM    PLC sensitivity to Ca2+
                                   "k_3k": 1.5e-03,  # /ms   IP3 Phosphorylation rate
                                   "k_5p": 0.01e-03,  # /ms   IP3 DePhosphorylation rate
                                   "Kf_plc": 0.35e-03,  # /ms   PLC-protein activation rate
                                   "Kb_plc": 2.2e-02,  # /ms   PLC-protein deactivation rate
                                   "Kf_gp": 0.33e-03,  # /ms   G-protein activation rate
                                   "Kb_gp": 2.17e-03,  # /ms   G-protein deactivation rate
                                   "delta_G": 0.01,  # G-protein intrinsic backgroung activity
                                   "Vr": 7.4,  # Maximal G-protein activation
                                   "Kr": 4467  # ug/mL A? concentration producing half-activation
                                   },
                            'ad': {"Vo": 0.19,  # uM    Intrinsic PLC-mediated IP3 production
                                        "Vq": 380,  # uM    Control parameter for influence of A? on IP3
                                        "Kip3k": 1.6,  # uM    Half-activation for 3-kinase
                                        "Kplc": 0.016,  # uM    PLC sensitivity to Ca2+
                                        "k_3k": 0.7e-03,  # /ms   IP3 Phosphorylation rate
                                        "k_5p": 0.005e-03,  # /ms   IP3 DePhosphorylation rate
                                        "Kf_plc": 0.75e-03,  # /ms   PLC-protein activation rate
                                        "Kb_plc": 2.0e-02,  # /ms   PLC-protein deactivation rate
                                        "Kf_gp": 0.047e-03,  # /ms   G-protein activation rate
                                        "Kb_gp": 4.7e-03,  # /ms   G-protein deactivation rate
                                        "delta_G": 0.012,  # G-protein intrinsic backgroung activity
                                        "Vr": 10,  # Maximal G-protein activation
                                        "Kr": 2000,  # ug/mL A? concentration producing half-activation}
                                        }
                            }

    release_params = {"allo": {"Kon": 0.097909,  # /uMms
                               "Koff": 3.316730,  # /ms
                               "I": 0.0000001,  # /ms
                               "F": 28.693830,
                               "b": 0.5008504,  # Cooperativity
                               "Krf": 1 / 6.339942,  # /ms rate of recovery of refractoriness
                               "Kmob": 0.003858,  # /uMms
                               "Kdemob": 0.002192,  # /ms
                               "Kprime": 0.028560,  # /uMms
                               "Kupr": 0.003124,  # /ms
                               "Kattach": 0.000144,  # /uMms
                               "Kdetach": 0.002413995},  # /ms
                      "dual": {
                          "alpha": 0.061200,  # /uMms   Association rate for synchronous release
                          "beta": 2.320000,  # /ms   Dissociation= rate for synchronous release
                          "Chi": 0.002933,  # /uMms   Association rate for Asynchronous release
                          "delta": 0.014829,  # /ms   Dissociation rate for Asynchronous release
                          "a": 0.025007,
                          "bo": 0.250007,  # Cooperative parameter
                          "gam_2": 2.000008,  # /ms Synchronous release rate
                          "gam_3": 0,  # /ms Asynchronous release rate
                          "gam_1": 0.000009,  # /ms Spontaneous release rate
                          "Krf2": 1 / 6.340000,  # /ms rate of recovery of refractoriness
                          "Kmob2": 0.000050,  # /uMms
                          "Kdemob2": 0.0022,  # /ms
                          "Kprime2": 0.027990,  # /uMms
                          "Kupr2": 0.005356,  # /ms
                          "Kattach2": 0.000150,  # /uMms
                          "Kdetach2": 0.001158}}  # /ms}

    release_params["dual"]["gam_3"] = release_params["dual"]["a"]*release_params["dual"]["gam_2"]

    for coupling_condition in coupling_conditions:


        data_directory = os.path.join("../../data/Train_data/" + coupling_condition + "/")

        if not os.path.exists(data_directory):
            os.mkdir(data_directory)

        for cell_condition in cell_conditions:

            for vgcc_channel_number in vgcc_channel_number_range:

                N_vgcc = vgcc_channel_number   # Number of vgccc

                DualSensorModel_release_rate = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1)))  # Overall release rate
                DualSensorModel_released_vesicles = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1)))  # Total vesicles released
                DualSensorModel_RRV = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Total Vesicles in the Release Ready Pool (i,e FRP and SRP)
                DualSensorModel_SpontRelRate = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Total Spontaneous vesicle release rate
                DualSensorModel_SyncRelRate = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) #  Total Synchronous vesicle release rate
                DualSensorModel_AsyncRelRate = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Total Aysnchronous vesicle release rate
                DualSensorModel_SlowRelRate = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Release rate of SRP
                DualSensorModel_FastRelRate = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Release rate of FRP
                DualSensorModel_Rate_ReservePool = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Rate of change of Reserve Pool (R)
                DualSensorModel_Rate_DockedPool = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Rate of change of Docked Pool (U)

                VGCC_Calcium = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Calcium concentration in AZ
                IP3R_Calcium = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Calcium concentration in IP3R cluster
                Cyto_Calcium = pd.DataFrame(index=list(range(1, num_trials+1)), columns=list(np.arange(1, nt+1))) # Calcium concentration in Cytosol


                for trial in range(1, num_trials+1):

                    [PoPQ, PoIPR, IcaPQ_type, time, y] = solve_ODEs(time_step, init, release_params, IP3Receptor_params,
                                                               volume_ratios, cell_condition, N_ipr, N_vgcc,IP3_params,
                                                               ISI, coupling_condition, dT_save, AP_type="paired pulse",
                                                               stimulus_strength=10)
                    VGCC_Calcium.loc[trial, :] = y[10, :]
                    IP3R_Calcium.loc[trial, :] = y[1, :]
                    Cyto_Calcium.loc[trial, :]= y[0, :]

### %%%%%%%%%%%%%%%% Allosteric Release Rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                    I = 0.0000001              # |     /ms     | vesicle fusion rate constant
#                   F = 28.693830              # |             | Vescile fusion cooperativity open Ca2+ binding

#                    AllostericModel_Slow.loc[trial, :] = I*(y[15,:] + F*y[16,:] + (F**2)*y[17,:] + (F**3)*y[18,:] + (F**4)*y[19,:] + (F**5)*y[20,:])
#                    AllostericModel_Fast.loc[trial, :] = I*(y[21,:] + F*y[22,:] + (F**2)*y[23,:] + (F**3)*y[24,:] + (F**4)*y[25,:] + (F**5)*y[26,:])
#                    AllostericModel.loc[trial, :] = AllostericModel_Slow.loc[trial, :] + AllostericModel_Fast.loc[trial, :]
#                    AllostericModel_release_probablity = scipy.integrate.cumtrapz(time, AllostericModel)
#                    AllostericModel_Cummulative = np.cumsum(AllostericModel)

### %%%%%%%%%%%%%%%% Dual Sensor Release Rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    gam_2 = release_params["dual"]["gam_2"]     # % /ms Synchronous release rate
                    gam_3 = release_params["dual"]["gam_3"]     # % /ms Asynchronous release rate
                    gam_1 = release_params["dual"]["gam_1"]     # % /ms Spontaneous release rate

                    R_dual, U_dual, RFv_dual, RFw_dual = y[27, :], y[28], y[29], y[30]

                    V00, V01, V02, V10, V11, V12, V20, V21, V22 = y[31, :], y[32, :], y[33, :], y[34, :], y[35, :], y[36, :], y[37, :], y[38, :], y[39, :]
                    V30, V31, V32, V40, V41, V42, V50, V51, V52 = y[40, :], y[41, :], y[42, :], y[43, :], y[44, :], y[45, :], y[46, :], y[47, :], y[48, :]

                    W00, W01, W02, W10, W11, W12, W20, W21, W22 = y[49, :], y[50, :], y[51, :], y[52, :], y[53, :], y[54, :], y[55, :], y[56, :], y[57, :]
                    W30, W31, W32, W40, W41, W42, W50, W51, W52 = y[58, :], y[59, :], y[60, :], y[61, :], y[62, :], y[63, :], y[64, :], y[65, :], y[66, :]

                    [alpha_dual, beta_dual, Chi_dual, delta_dual, a_dual, bo_dual, gam_2_dual, gam_3_dual, gam_1_dual,
                     Krf_dual, Kmob_dual, Kdemob_dual, Kprime_dual, Kupr_dual, Kattach_dual, Kdetach_dual] = release_params['dual'].values()

                    DualSensorModel_Rate_ReservePool.loc[trial, :] = -Kmob_dual*Cyto_Calcium.loc[trial, :]*R_dual + U_dual*Kdemob_dual
                    DualSensorModel_Rate_DockedPool.loc[trial, :] = -Kprime_dual*U_dual*Cyto_Calcium.loc[trial, :]*(1 - RFv_dual) + Kupr_dual*V00


                    DualSensorModel_SlowRelRate.loc[trial,:] = gam_3*(V02 + V12 + V22 + V32 + V42 + V52) + gam_2*(V50 + V51 + V52) + gam_1*V00
                    DualSensorModel_FastRelRate.loc[trial,:] = gam_3*(W02 + W12 + W22 + W32 + W42 + W52) + gam_2*(W50 + W51 + W52) + gam_1*W00

                    DualSensorModel_release_rate.loc[trial, :] = DualSensorModel_SlowRelRate.loc[trial,:] + \
                                                                 DualSensorModel_FastRelRate.loc[trial,:]

                    DualSensorModel_released_vesicles.loc[trial, :] = scipy.integrate.cumtrapz(time, DualSensorModel_release_rate.loc[trial, :])
                    num_vesicles_in_FRP = 5
                    num_vesicles_in_SRP = 5

                    DualSensorModel_Slow_RRV = V00 + V01 + V02 + V10 + V11 + V12 + V20 + V21 + V22 + \
                                                     V30 + V31 + V32 + V40 + V41 + V42 + V50 + V51 + V52
                    DualSensorModel_Fast_RRV =  W00 + W01 + W02 + W10 + W11 + W12 + W20 + W21 + W22 + \
                                                     W30 + W31 + W32 + W40 + W41 + W42 + W50 + W51 + W52


                    DualSensorModel_RRV.loc[trial, :] = DualSensorModel_Slow_RRV + DualSensorModel_Fast_RRV


                    DualSensorModel_AsyncRelRate.loc[trial, :] = gam_3*(V02 + V12 + V22 + V32 + V42 + V52) + \
                                                                gam_3*(W02 + W12 + W22 + W32 + W42 + W52)
                    DualSensorModel_SpontRelRate.loc[trial, :] = gam_1*(V00 + W00)
                    DualSensorModel_SyncRelRate.loc[trial, :] = gam_2*(V50 + V51 + V52) + gam_2*(W50 + W51 + W52)


                filename  = f"{data_directory}/{cell_condition}_Train_Rate_ReservePool_{N_vgcc}_VGCC.txt"
                DualSensorModel_Rate_ReservePool.mean(axis=0).to_csv(filename, header=False, index=False)

                filename =  f"{data_directory}/{cell_condition}_Train_Rate_DockedPool_{N_vgcc}_VGCC.txt"
                DualSensorModel_Rate_DockedPool.mean(axis=0).to_csv(filename, header=False, index=False)

                filename = f"{data_directory}/{cell_condition}_Train_AsyncRelRate_{N_vgcc}_VGCC.txt"
                DualSensorModel_AsyncRelRate.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  = f"{data_directory}/{cell_condition}_Train_SyncRelRate_{N_vgcc}_VGCC.txt"
                DualSensorModel_SyncRelRate.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  = f"{data_directory}/{cell_condition}_Train_SpontRelRate_{N_vgcc}_VGCC.txt"
                DualSensorModel_SpontRelRate.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  = f"{data_directory}/{cell_condition}_Train_SlowRelRate_{N_vgcc}_VGCC.txt"
                DualSensorModel_SlowRelRate.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  =  f"{data_directory}/{cell_condition}_Train_FastRelRate_{N_vgcc}_VGCC.txt"
                DualSensorModel_FastRelRate.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  = f"{data_directory}/{cell_condition}_Train_RelVes_{N_vgcc}_VGCC.txt"
                DualSensorModel_released_vesicles.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  = f"{data_directory}/{cell_condition}_Train_RelRate_{N_vgcc}_VGCC.txt"
                DualSensorModel_release_rate.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  = f"{data_directory}/{cell_condition}_Train_RRV_{N_vgcc}_VGCC.txt"
                DualSensorModel_RRV.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  = f"{data_directory}/{cell_condition}_Train_IP3_Calcium_{N_vgcc}_VGCC.txt"
                IP3R_Calcium.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  = f"{data_directory}/{cell_condition}_Train_VGCC_Calcium_{N_vgcc}_VGCC.txt"
                VGCC_Calcium.mean(axis=0).to_csv(filename, header=False, index=False)

                filename  = f"{data_directory}/{cell_condition}_Train_CYTO_Calcium_{N_vgcc}_VGCC.txt"
                Cyto_Calcium.mean(axis=0).to_csv(filename, header=False, index=False)




