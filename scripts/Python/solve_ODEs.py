import numpy as np
from kinetic_schemes import kinetic_schemes
from ODEs import ODEs


def solve_ODEs(time_step, init, release_params, IP3Receptor_params,
               volume_ratios, cell_condition, N_ipr, N_vgcc,IP3_params,
               ISI, coupling_condition, dT_save, AP_type,
               stimulus_strength):

    dt = time_step
    nn = int(dT_save / dt)  # skip every nn steps before saving data

    # # Initial Conditions For Stochastic Simulation
    state_ip3r = np.ones((1, N_ipr))  # initial state of each channel. R = 1, O = 2, I = 3
    state_PQ = np.zeros((1, N_vgcc))  # Initial state of each VGCC. C = 0, O = 1

    num_ode = len(init)  # number of odes
    y = np.zeros((num_ode, time_step))  # save all avriables for plotting

    for i in range(len(init)):
        y[i, 0] = init[i]

    IcaPQ_type = np.zeros((1, time_step))
    PoIPR = np.zeros((1, time_step))
    PoPQ = np.zeros((1, time_step))
    tt = 0  # initialize time t = 0
    t = np.zeros((1, time_step)) # save time

    for i in range(1, time_step):
        yn = y[:, i - 1]  # variable values at previous time step

        # Runge-Kutta Method
        for j in range(0, nn):  # this loop is used so that we don't have to save data at each time step

            [Po_PQ, Po_ip3r, IcaPQ, state_ip3r, state_PQ] = kinetic_schemes(tt, yn, dt, state_ip3r, state_PQ,
                                                                           N_vgcc, N_ipr, IP3Receptor_params,
                                                                           cell_condition
                                                                           )

            k1 = ODEs(tt, yn,  Po_ip3r, IcaPQ, volume_ratios, AP_type, stimulus_strength, ISI, N_vgcc,
                        coupling_condition, cell_condition, IP3_params, release_params
                        )
            k2 = ODEs(tt, yn + [k1i * (dt / 2) for k1i in k1],   Po_ip3r, IcaPQ, volume_ratios, AP_type, stimulus_strength,
                        ISI, N_vgcc, coupling_condition, cell_condition, IP3_params, release_params
                        )
            k3 = ODEs(tt, yn + [k2i * (dt / 2) for k2i in k2], Po_ip3r, IcaPQ, volume_ratios, AP_type, stimulus_strength,
                        ISI, N_vgcc, coupling_condition, cell_condition, IP3_params, release_params
                        )
            k4 = ODEs(tt, yn + [k3i * dt for k3i in k3],   Po_ip3r, IcaPQ, volume_ratios, AP_type, stimulus_strength,
                        ISI, N_vgcc, coupling_condition, cell_condition, IP3_params, release_params
                      )
            k_in = [k1[i] + [2 * k2i for k2i in k2][i] + [2 * k3i for k3i in k3][i] + k4[i] for i in range(len(k4))]
            yn = yn + [k_ini * (dt / 6) for k_ini in k_in]

            PoPQ[i, 0] = Po_PQ
            PoIPR[i, 0] = Po_ip3r
            IcaPQ_type [i, 0] = IcaPQ

        y[:, i] = yn  # update the new value of all variables
        tt = i * dT_save  # generate time array for plotting
        t[i] = tt
    return [PoPQ, PoIPR, IcaPQ_type, t, y]
