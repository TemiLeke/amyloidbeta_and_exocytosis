clear all; close all;

%% Description

%       This script simulates calcium dynamics, neurotransmitter release in
%       physiologically plausible hippocampal CA1 terminal model. 
%       All simulations have a time step of 1 us (i,e 0.001 ms). We
%       simulated fo 50 trials and calculated average response for VGGC 
%       values in the range= 5--150 (represented as k in simulation). 
%       Release probability in response to an AP is calculated by
%       integrating the release rate for slow and fast 
%       vesicles divided by the number of RRV initially in both pools.
%       Simulations for single AP and paired pulse protocol are run for
%       100ms in seperate scripts, while AP Train simulation were done for 450 ms. 
%       For paired pulse protocol, the total duration of each pulse is 30ms. ISI can be
%       modified as desired. AP train of 20Hz was used in our simulations.
%       Simulation results (average of 50 trials) are saved to an automatically created data
%       directory for each VGCC number.
 

%% Model Parameters and Configuration

%       All kinetic and reaction rates/parameters are specified in the
%       kinetic_schemes.m and ODEs.m file. Each paramter is annotated
%       accordingly. 
%       All state variables are set to closed, and all dynamical variables
%       are set to physiologically plausible values. Current stimulation is
%       only applied 3 ms after initialization of simulation. This ensures
%       that all varaibles attain thier steady-state resting levels. 
%       The coupled ODE system is solved using the RK4 algorithm as other
%       inbuilt methods appear unstable with stochastic kinetic schemes.


%% Note:

%       Please note that all files; "kinetic_schemes.m", "ODEs.m",
%       "solve_Train_ODEs.m", "solve_PPR_ODEs.m", and
%       "solveSingleAP_ODEs.m" have to be in the same directory for
%       simulation to run without hiccups. 


%% Simulation

global dt N_ip3r N_vgcc N_vesicles N_pores state_ip3r state_PQ time SCL; % set global variables


coupling_condition = "Same_Coupling";               %  Coupling configuration
ca_source = ["abeta", "ip3r_nc", "ip3r_nc_and_abeta", "ip3r_hc_and_abeta"];    %  Additional Ca Source to the AZ

for index=1:length(ca_source)
    calcium_source = ca_source(index);

    if (calcium_source == "ip3r_hc_and_abeta" || calcium_source=="abeta")
        cell_condition = "AD";                                                  %  Cell condition (WT or AD)
    elseif (calcium_source == "ip3r_nc_and_abeta" || calcium_source=="ip3r_nc") 
        cell_condition = "WT";                                                  %  Cell condition (WT or AD) 
    end



    data_directory = strcat("../../data/APTrain_data_mop/", calcium_source, "/");
    mkdir(data_directory)

    MOP = 0.015;                      %    Mean open probability of the pore at 20 minutes
    n_mot_percent_incr = 100;         %    average percentage increase in mean open probability (MOP).        
    
    for time_of_recording=20:20:80      % time of recording in minutes.      

        for k=35     % Interation over Number of channels ranging from 5 to 150 VGCC number

            num_trials = 50;                                                        % Number of trials
            time_steps = 450000;           % numbser of time steps     (equivalent to 100 ms for steps of 0.001 ms)
            % Matrix for storing release model variables for all trials             % Dimension --> (Number of trails, Number of time steps)

            DualSensorModel_release_rate = zeros(num_trials, time_steps);           % Overall release rate 
            DualSensorModel_SpontRelRate = zeros(num_trials, time_steps);          % Total Spontaneous vesicle release rate  
            DualSensorModel_SyncRelRate = zeros(num_trials, time_steps);           % Total Synchronous vesicle release rate  
            DualSensorModel_AsyncRelRate = zeros(num_trials, time_steps);          % Total Aysnchronous vesicle release rate  
            DualSensorModel_SlowRelRate = zeros(num_trials, time_steps);           % Release rate of SRP  
            DualSensorModel_FastRelRate = zeros(num_trials, time_steps);           % Release rate of FRP  
            DualSensorModel_Slow_RRV = zeros(num_trials, time_steps);              % Total Vesicles in the Slow Release Ready Pool (i,e SRP)
            DualSensorModel_Fast_RRV = zeros(num_trials, time_steps);              % Total Vesicles in the Fast Release Ready Pool (i,e FRP)
            ActionPotential_Train = zeros(num_trials, time_steps);         

            DualSensorModel_released_vesicles = zeros(num_trials, time_steps);     % Total vesicles released 
            DualSensorModel_RRV = zeros(num_trials, time_steps);                   % Total Vesicles in the Release Ready Pool (i,e FRP and SRP)    
            DualSensorModel_Rate_ReservePool = zeros(num_trials, time_steps);      % Rate of change of Reserve Pool (R)  
            DualSensorModel_Rate_DockedPool = zeros(num_trials, time_steps);       % Rate of change of Docked Pool (U)

            VGCC_Calcium = zeros(num_trials, time_steps);       % Calcium concentration in AZ
            IP3R_Calcium = zeros(num_trials, time_steps);       % Calcium concentration in IP3R cluster
            Cyto_Calcium = zeros(num_trials, time_steps);       % Calcium concentration in Cytosol    

            for n=1:num_trials
                %% Parameters for Simulation

                nt = time_steps;       % numbser of time steps     (equivalent to 100 ms for steps of 0.001 ms)
                dt = 0.001;            % time step for integration (1 uS (micro-second) or 0.001 mS) 
                dT_save = 0.001;       % time step for saving data
                nn = fix(dT_save/dt);  % skip every nn steps before saving data
                time = 0;              % initial time            
                num_ode = 67;          % Total number of all odes

                %% Initial Conditions For Stochastic Simulation

                N_ip3r  = 10;                        % number of IP3R channels 
                N_vgcc = k;                          % number of VGCC (k is specified in loop)
                N_vesicles = 200;                    % number of vesicles
                N_pores = 200;                       % number of abeta pores
                SCL = 1;                             % Subconductance level
                state_ip3r = ones(1, N_ip3r);        % initial state of each IP3R channel is R. Markov Chain Legend R=1, A=2, O=3, I=4
                state_PQ = zeros(1, N_vgcc);         % Initial state of each VGCC ia C1. Markov Chain Legend C1=1, C2=2, C3=3, C4=4, O=5

                %% Initial Conditions For Deterministic Simulation

                y = zeros(num_ode, nt);              % Solutions. dimension = (Number of ODes x Number of Time Steps)

                y(:,1) = [0.1; 0.1; 56; -70; 0.01; 0.01;...  % dCacdt; dCa_iprdt; dCtdt; dvdt; dn2dt; dh2dt
                          1; 1; 0.16; 0.1;...                % dPLCdt; dGdt; dIP3dt; dCabetadt;
                          0.1;...                            % dCa_vgccdt
                          170; 20; 0; 0; 5; 0; 0; 0; 0; 0; 5; 0; 0; 0; 0; 0;...    % dR_allodt; dU_allodt; dRFv_allodt; dRFw_allodt; dV0dt; dV1dt; dV2dt; dV3dt; dV4dt; dV5dt; dW0dt; dW1dt; dW2dt; dW3dt; dW4dt; dW5dt
                          170; 20; 0; 0;...                                        % dR_dualdt; dU_dualdt; dRFv_dualdt; dRFw_dualdt
                          5; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;... % dV00; dV01; dV02; dV10; dV11; dV12; dV20; dV21; dV22; dV30; dV31; dV32; dV40; dV41; dV42; dV50; dV51; dV52;
                          5; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];  % dW00; dW01; dW02; dW10; dW11; dW12; dW20; dW21; dW22; dW30; dW31; dW32; dW40; dW41; dW42; dW50; dW51; dW52;


                %% Simulation

                AP_condition = "Train";
                ISI = 40;                                           % ms    Inter-spike Interval for Paired-Pulse protocol
                abeta_dose = 0;                                     % ug/mL Amyloid beta dose used WT: 0 ug/mL AD: 3 ug/mL
                IcaPQ_type = zeros(1, nt);                          % uA    Initialize all Single channel currents to zero 
                Pr_ip3r = zeros(1, nt);                             %       Initialize all IP3R open probabilities to zero 
                t = 0;                                              % ms    Initialize time t = 0
                
                init_open_time = 0;                                 % ms    Initial abeta open time
                close_time = nt*dt;                                 % ms    Time when all pores close
                mean_open_time = 10;                                % ms    Mean open duration of the pore
                exp_duration = nt*dt;                               % ms    duration of experimental recordings corresponding to (20-30s) in Syed Paper ; https://doi.org/10.1101/2022.04.29.490101 doi:               
                mean_open_proba = MOP;                              %       Mean open probability at time of recording
                total_open_time = mean_open_proba*exp_duration;     % ms    Total duration of pore opening 

                if  total_open_time/mean_open_time < 1
                    sample_times = init_open_time:dt:close_time;
                    num_open_events = 1;                                                    %   number of events
                    open_times = sort(randsample(sample_times, num_open_events));           %   randomly chosen open times of the abeta pores
                    pore_open_times = (open_times(1):dt:(open_times(1)+total_open_time));
                else 
                    times = init_open_time:dt:close_time;
                    sample_times = times(1:mean_open_time/dt:end);
                    num_open_events = round(total_open_time/mean_open_time);                       %   number of events
                    open_times = sort(randsample(sample_times, num_open_events));           %   randomly chosen open times of the abeta pores
                    pore_open_times = (open_times(1):dt:(open_times(1)+mean_open_time));

                    for open_event=2:length(open_times)
                        pore_open_times = [pore_open_times, (open_times(open_event):dt:(open_times(open_event) + mean_open_time))];
                    end

                end
                

                for i = 2:nt
                    yn = y(:,i-1);                                         % variable values at previous time step  

                    % Runge-Kutta Method

                    for j = 1:nn                                                     % this loop is used so that we don't have to save data at each time step
                        [Po_ip3r, IcaPQ] = kinetic_schemes(yn, cell_condition);      % compute the open probability at current time   
                        k1 = ODEs(t, yn, Po_ip3r, IcaPQ, AP_condition, cell_condition, ISI, abeta_dose, coupling_condition, pore_open_times, calcium_source);                          
                        k2 = ODEs(t, yn + k1*dt/2, Po_ip3r, IcaPQ, AP_condition, cell_condition, ISI, abeta_dose, coupling_condition, pore_open_times, calcium_source);
                        k3 = ODEs(t, yn + k2*dt/2, Po_ip3r, IcaPQ, AP_condition, cell_condition, ISI, abeta_dose, coupling_condition, pore_open_times, calcium_source);
                        k4 = ODEs(t, yn + k3*dt, Po_ip3r, IcaPQ, AP_condition, cell_condition, ISI, abeta_dose, coupling_condition, pore_open_times, calcium_source);
                        yn = yn + (k1 + 2*k2+ 2*k3 + k4)* dt/6;                      % Obtain solutions to ODEs at current time 
                    end

                    y(:, i) = yn;                   % update the new value of all variables
                    Pr_ip3r(i) = Po_ip3r;           % save the open probability at current time
                    IcaPQ_type(i) = IcaPQ;          % save the single channel current at current time
                    t = (i)*dT_save;                % generate time array

                end

                t = (1: nt)*dT_save;        % generate time array for plotting


                VGCC_Calcium(n, :) = y(11, :);
                IP3R_Calcium(n, :) = y(2, :);
                Cyto_Calcium(n, :) = y(1, :);

                %% Allosteric Release Rates

                I = 0.0000001;              % |     /ms     | vesicle fusion rate constant
                F = 28.693830;              % |             | Vescile fusion cooperativity open Ca2+ binding

                AllostericModel_Slow(1,:) = I.*(y(16,:) + F.*y(17,:) + (F^2).*y(18,:) + (F^3).*y(19,:) + (F^4).*y(20,:) + (F^5).*y(21,:));          
                AllostericModel_Fast(1,:) = I.*(y(22,:) + F.*y(23,:) + (F^2).*y(24,:) + (F^3).*y(25,:) + (F^4).*y(26,:) + (F^5).*y(27,:));     
                AllostericModel(1,:) = AllostericModel_Slow(1,:) + AllostericModel_Fast(1,:);
                AllostericModel_release_probablity = cumtrapz(t, AllostericModel);
                AllostericModel_Cummulative = cumsum(AllostericModel);


                %% Dual Sensor Release Rates

                a = 0.025007;
                gam_2 = 2.000008;     % /ms Synchronous release rate
                gam_3 = a*gam_2;      % /ms Asynchronous release rate
                gam_1 = 0.000009;     % /ms Spontaneous release rate

                R2 = y(28, :) ; U2 = y(29, :) ; RFv2 = y(30, :) ; RFw2 = y(31, :) ;
                V00 = y(32,:); V01 = y(33,:); V02 = y(34,:); V10 = y(35,:); V11 = y(36,:); V12 = y(37,:); V20 = y(38,:); V21 = y(39,:); V22 = y(40,:);... 
                V30 = y(41,:); V31 = y(42,:); V32 = y(43,:); V40 = y(44,:); V41 = y(45,:); V42 = y(46,:); V50 = y(47,:); V51 = y(48,:); V52 = y(49,:);

                W00 = y(50,:); W01 = y(51,:); W02 = y(52,:); W10 = y(53,:); W11 = y(54,:); W12 = y(55,:); W20 = y(56,:); W21 = y(57,:); W22 = y(58,:); 
                W30 = y(59,:); W31 = y(60,:); W32 = y(61,:); W40 = y(62,:); W41 = y(63,:); W42 = y(64,:); W50 = y(65,:); W51 = y(66,:); W52 = y(67,:);




            % Reccruitment model parameters |   Unit    |   Description
                Krf_dual =  1/10.340000;  % |   /ms     |   rate of recovery of refractoriness
                Kmob_dual = 0.000050;     % |   /uMms   |   Mobilization rate
                Kdemob_dual = 0.0022;     % |   /ms     |   Demobilization rate  
                Kprime_dual = 0.027990;   % |   /uMms   |   Priming rate 
                Kupr_dual = 0.005356;     % |   /ms     |   Unpriming rate
                Kattach_dual = 0.00150;   % |   /uMms   |   Attachement rate
                Kdetach_dual = 0.001158;  % |   /ms     |   Detachhement rate


                DualSensorModel_ReservePool(1, :) = -Kmob_dual.*Cyto_Calcium(n, :).*R2 + U2.*Kdemob_dual;
                DualSensorModel_DockedPool(1, :) = -Kprime_dual.*U2.*Cyto_Calcium(n, :).*(1 - RFv2) + Kupr_dual.*V00;
                DualSensorModel_Rate_ReservePool(n, :) = DualSensorModel_ReservePool(1, :);
                DualSensorModel_Rate_DockedPool(n, :) = DualSensorModel_DockedPool(1, :);

                DualSensorModel_Slow(1,:) = gam_3.*(V02 + V12 + V22 + V32 + V42 + V52) + ...
                                            gam_2.*(V50 + V51 + V52) + gam_1.*V00;

                DualSensorModel_Fast(1,:) = gam_3.*(W02 + W12 + W22 + W32 + W42 + W52) + ...
                                            gam_2.*(W50 + W51 + W52) + gam_1.*W00; 

                DualSensorModel(1,:) = DualSensorModel_Slow + DualSensorModel_Fast;

                DualSensorModel_release_rate(n, :) = DualSensorModel(1, :);
                DualSensorModel_released_vesicles(n,:) = cumtrapz(t, DualSensorModel(1, :));

                DualSensorModel_FastRelRate(n, :) = DualSensorModel_Fast(1,:);
                DualSensorModel_SlowRelRate(n, :) = DualSensorModel_Slow(1,:);


                DualSensorModel_SlowRRV(1,:) = V00 + V01 + V02 + V10 + V11 + V12 + V20 + V21 + V22 + ...
                                                V30 + V31 + V32 + V40 + V41 + V42 + V50 + V51 + V52;

                DualSensorModel_FastRRV(1,:) =  W00 + W01 + W02 + W10 + W11 + W12 + W20 + W21 + W22 + ...
                                                 W30 + W31 + W32 + W40 + W41 + W42 + W50 + W51 + W52;

                DualSensorModelRRV(1, :) = DualSensorModel_SlowRRV + DualSensorModel_FastRRV;

                DualSensorModel_RRV(n,:) = DualSensorModelRRV(1, :);

                DualSensorModel_Slow_RRV(n,:) = DualSensorModel_SlowRRV(1, :);
                DualSensorModel_Fast_RRV(n,:) = DualSensorModel_FastRRV(1, :);

                DualSensor_Spontaneous(1,:) = gam_1.*(V00 + W00);
                DualSensor_Synchronous(1,:) =  gam_2.*(V50 + V51 + V52) + gam_2*(W50 + W51 + W52); 
                DualSensor_Asynchronous(1,:) = gam_3.*(V02 + V12 + V22 + V32 + V42 + V52) + ...
                                               gam_3.*(W02 + W12 + W22 + W32 + W42 + W52);


                DualSensorModel_AsyncRelRate(n, :) = DualSensor_Asynchronous(1,:);  
                DualSensorModel_SpontRelRate(n, :) = DualSensor_Spontaneous(1,:); 
                DualSensorModel_SyncRelRate(n, :) = DualSensor_Synchronous(1,:); 


            end

        %% Saving Data


            filename  = sprintf("Train_Rate_ReservePool_%d.txt", k);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_Rate_ReservePool, 1));
            fclose(fileID);

            filename  = sprintf("Train_Rate_DockedPool_%d.txt", k);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_Rate_DockedPool, 1));
            fclose(fileID);

            filename  = sprintf("Train_AsyncRelRate_%d.txt", k);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_AsyncRelRate, 1));
            fclose(fileID);

            filename  = sprintf("Train_SyncRelRate_%d.txt", k);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_SyncRelRate, 1));
            fclose(fileID);

            filename  = sprintf("Train_SpontRelRate_%d.txt", k);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_SpontRelRate, 1));
            fclose(fileID);

            filename  = sprintf("Train_SlowRelRate_%d.txt", k);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_SlowRelRate, 1));
            fclose(fileID);

            filename  = sprintf("Train_FastRelRate_%d.txt", k);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_FastRelRate, 1));
            fclose(fileID);

            filename  = sprintf("Train_RelRate_%d_min.txt", time_of_recording);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_release_rate, 1));
            fclose(fileID);

            filename  = sprintf("Train_RelVes_%d_min.txt", time_of_recording);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_released_vesicles, 1));
            fclose(fileID);

            filename  = sprintf("Train_RRV_%d_min.txt", time_of_recording);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_RRV, 1));
            fclose(fileID);  

            filename  = sprintf("Train_Fast_RRV_%d_min.txt", time_of_recording);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_Fast_RRV, 1));
            fclose(fileID); 

            filename  = sprintf("Train_Slow_RRV_%d_min.txt", time_of_recording);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(DualSensorModel_Slow_RRV, 1));
            fclose(fileID); 

            filename  = sprintf("Train_IP3_Calcium_%d_min.txt", time_of_recording);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(IP3R_Calcium, 1));
            fclose(fileID);

            filename  = sprintf("Train_VGCC_Calcium_%d_min.txt", time_of_recording);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(VGCC_Calcium, 1));
            fclose(fileID);

            filename  = sprintf("Train_CYTO_Calcium_%d_min.txt", time_of_recording);
            fileID = fopen(strcat(data_directory, filename), "w");
            fprintf(fileID, "%.15f\n", mean(Cyto_Calcium, 1));
            fclose(fileID);

            
        end
       
        MOP = MOP + MOP;                                    %       save MOP for next iteration/time of recording

    end
end