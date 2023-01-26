 clear all; close all; 


%% Data and figure  directory

mkdir ../../results Figure_1

root = strcat(fileparts(fileparts(pwd)), "/data/SingleAP_data/");

%%   PPR, Probability, Release Rates and other plots

base_dir = strcat(fileparts(fileparts(pwd)), "/data/");
t = importdata(strcat(base_dir, "singleap_time.txt"));

num_conditions = 3;                      % condition 1 corresponds to SamCoupling, 2 - AD Higher Coupling, 3 - WT Higher Coupling

total_rel_ves_WT = cell(1, num_conditions);
total_RRV_WT = cell(1, num_conditions);
rel_proba_WT = cell(1, num_conditions);
rel_rate_WT = cell(1, num_conditions);
Sync_RelRate_WT = cell(1, num_conditions);
Async_RelRate_WT = cell(1, num_conditions);
Spont_RelRate_WT = cell(1, num_conditions);
Docked_Rate_WT = cell(1, num_conditions);
Reserved_Rate_WT = cell(1, num_conditions);
Slow_Rate_WT = cell(1, num_conditions);
Fast_Rate_WT = cell(1, num_conditions);
Ca_VGCC_WT = cell(1, num_conditions);
Ca_CYTO_WT = cell(1, num_conditions);
Ca_IP3R_WT = cell(1, num_conditions); 
Ca_ABETA_WT = cell(1, num_conditions); 

min_channel_num = 5;
max_channel_num = 150;
num_channels = int32(((max_channel_num - 5)/5) + 1);
coupling_conditions = ["ip3r_nc", "ip3r_nc_and_abeta", "ip3r_hc_and_abeta"];


for coupling_cond=1:length(coupling_conditions)


    condition = coupling_conditions(coupling_cond);
    

    if condition == "ip3r_nc"
        k = 1;
    elseif condition == "ip3r_nc_and_abeta"
        k = 2;
    elseif condition == "ip3r_hc_and_abeta"
        k = 3;
    end


    total_rel_ves_WT{k} = zeros(length(t), num_channels);
    total_RRV_WT{k} = zeros(length(t), num_channels);
    rel_proba_WT{k} = zeros(length(t), num_channels);
    rel_rate_WT{k} = zeros(length(t), num_channels);
    Sync_RelRate_WT{k} = zeros(length(t), num_channels);
    Async_RelRate_WT{k} = zeros(length(t), num_channels);
    Spont_RelRate_WT{k} = zeros(length(t), num_channels);
    Docked_Rate_WT{k} = zeros(length(t), num_channels);
    Reserved_Rate_WT{k} = zeros(length(t), num_channels);
    Slow_Rate_WT{k} = zeros(length(t), num_channels);
    Fast_Rate_WT{k} = zeros(length(t), num_channels);
    Ca_VGCC_WT{k} = zeros(length(t), num_channels);
    Ca_CYTO_WT{k} = zeros(length(t), num_channels);
    Ca_IP3R_WT{k} = zeros(length(t), num_channels);
    Ca_ABETA_WT{k} = zeros(length(t), num_channels);
    

    for channel_number=5:5:150
        j = int32(((channel_number - 5)/5) + 1);

        total_rel_ves_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_RelVes_", num2str(channel_number), "_VGCC", ".txt"));
        total_RRV_WT{k}(:, j) = importdata(strcat(root, condition, '/',"SingleAP_RRV_", num2str(channel_number), "_VGCC",  ".txt"));
        rel_proba_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_RelProba_", num2str(channel_number), "_VGCC",  ".txt"));
        rel_rate_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_RelRate_", num2str(channel_number), "_VGCC", ".txt"));
        Sync_RelRate_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_SyncRelRate_", num2str(channel_number), "_VGCC", ".txt"));
        Async_RelRate_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_AsyncRelRate_", num2str(channel_number), "_VGCC", ".txt"));
        Spont_RelRate_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_SpontRelRate_", num2str(channel_number), "_VGCC", ".txt"));
        Docked_Rate_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_Rate_DockedPool_", num2str(channel_number), "_VGCC", ".txt"));
        Reserved_Rate_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_Rate_ReservePool_", num2str(channel_number), "_VGCC", ".txt"));
        Slow_Rate_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_SlowRelRate_", num2str(channel_number), "_VGCC", ".txt"));
        Fast_Rate_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_FastRelRate_", num2str(channel_number), "_VGCC", ".txt"));
        Ca_VGCC_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_VGCC_Calcium_", num2str(channel_number), "_VGCC", ".txt"));
        Ca_CYTO_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_CYTO_Calcium_", num2str(channel_number), "_VGCC", ".txt"));
        Ca_IP3R_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_IP3_Calcium_", num2str(channel_number), "_VGCC", ".txt")); 
        Ca_ABETA_WT{k}(:, j) = importdata(strcat(root, condition, '/', "SingleAP_ABETA_Calcium_", num2str(channel_number), "_VGCC", ".txt")); 

  
    end
end

%% params

num_channels = int32(((150 - 5)/5) + 1);
channel_number = min_channel_num:5:max_channel_num;

%% Release Probability versus channel number For all calcium source combination. 

Pr_matrix_WT = cell(1, num_conditions);

params_init_pr = [0.005 1.127 1.903 2.716];
lb_pr = [];
ub_pr = [];

params_pr_WT = cell(1, num_conditions);
Pr_matrix_WT_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end

    %Pr_matrix_WT{condition_index} = rel_proba_WT{condition_index}(end, :);

    Pr_matrix_WT{condition_index} = total_rel_ves_WT{condition_index}(30000, :)./total_RRV_WT{condition_index}(30000-1, :);
    
    x = log10(channel_number);
    Y_WT = Pr_matrix_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)DoseResponseFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);
                                    
    params_pr_WT{condition_index} = params_WT;
    
    Pr_matrix_WT_fit{condition_index} = DoseResponseFit(params_pr_WT{condition_index}, log10(channel_number));
    
end


%% Peak Release Rate versus channel number

Peak_RelRate_WT = cell(1, num_conditions);
Peak_RelRate_AD = cell(1, num_conditions);

params_init_prr = [6.848 0.03423 1.974 -2.547];
lb_prr = [];
ub_prr = [];

params_prr_WT = cell(1, num_conditions);
params_prr_AD = cell(1, num_conditions);
Peak_RelRate_WT_fit = cell(1, num_conditions);
Peak_RelRate_AD_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end

        
    Peak_RelRate_WT{condition_index} = max(rel_rate_WT{condition_index}, [], 1);    
    
    x = log10(channel_number);
    Y_WT = Peak_RelRate_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)DoseResponseFit(params_WT, x),params_init_prr, x, Y_WT, lb_prr,ub_prr);
                                    
    params_prr_WT{condition_index} = params_WT;
    
    Peak_RelRate_WT_fit{condition_index} = DoseResponseFit(params_prr_WT{condition_index}, log10(channel_number));
   
    
end


%% Peak Synchronous Release Rate versus channel number

SyncPeak_RelRate_WT = cell(1, num_conditions);
SyncPeak_RelRate_AD = cell(1, num_conditions);

params_init_prr = [0.03352 6.844 1.974  2.547];
lb_prr = [];
ub_prr = [];

params_prr_WT = cell(1, num_conditions);
params_prr_AD = cell(1, num_conditions);
SyncPeak_RelRate_WT_fit = cell(1, num_conditions);
SyncPeak_RelRate_AD_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end
 
    SyncPeak_RelRate_WT{condition_index} = max(Sync_RelRate_WT{condition_index}, [], 1);
    
    
    x = log10(channel_number);
    Y_WT = SyncPeak_RelRate_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)DoseResponseFit(params_WT, x),params_init_prr, x, Y_WT, lb_prr,ub_prr);
                                    
    params_prr_WT{condition_index} = params_WT;
    
    SyncPeak_RelRate_WT_fit{condition_index} = DoseResponseFit(params_prr_WT{condition_index}, log10(channel_number));
   
    
end


%% Peak Asynchronous Release Rate versus channel number

AsyncPeak_RelRate_WT = cell(1, num_conditions);
AsyncPeak_RelRate_AD = cell(1, num_conditions);

params_init_prr = [3.32e-08 5.058  0.001488 0.1431];
lb_prr = [];
ub_prr = [];

params_prr_WT = cell(1, num_conditions);
params_prr_AD = cell(1, num_conditions);
AsyncPeak_RelRate_WT_fit = cell(1, num_conditions);
AsyncPeak_RelRate_AD_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end

        
    AsyncPeak_RelRate_WT{condition_index} = max(Async_RelRate_WT{condition_index}, [], 1);
    
    
    x = log10(channel_number);
    Y_WT = AsyncPeak_RelRate_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_prr, x, Y_WT, lb_prr,ub_prr);
                                    
    params_prr_WT{condition_index} = params_WT;
    
    AsyncPeak_RelRate_WT_fit{condition_index} = BiexponentialFit(params_prr_WT{condition_index}, log10(channel_number));
   
    
end

%% Cumulative AZ Calcium Concentration (100 ms) after stimulation Versus Release Probability 

Cummulative_VGCC_Ca_WT = cell(1, num_conditions);
Cummulative_VGCC_Ca_AD = cell(1, num_conditions);


params_init_cum = [18.21 6.911 0.2096 ];
lb_cum = [];
ub_cum = [];

params_cum_WT = cell(1, num_conditions);
params_cum_AD = cell(1, num_conditions);
Cummulative_VGCC_Ca_WT_fit = cell(1, num_conditions);
Cummulative_VGCC_Ca_AD_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end


    VGCC_Ca_WT = cumtrapz(t, Ca_VGCC_WT{condition_index}, 1);
     
    Cummulative_VGCC_Ca_WT{condition_index} = VGCC_Ca_WT(100000,:);
    
    x_WT = Pr_matrix_WT{condition_index};
    Y_WT = Cummulative_VGCC_Ca_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)LogitFit(params_WT, x_WT),params_init_cum, x_WT, Y_WT, lb_cum,ub_cum);
                                    
    params_cum_WT{condition_index} = params_WT;
    
    Cummulative_VGCC_Ca_WT_fit{condition_index} = LogitFit(params_cum_WT{condition_index},Pr_matrix_WT{condition_index});
    
       
end

%% Total Vesicle Released (100 ms after stimulation) vs Release Probability


Vesicles_Released_WT = cell(1, num_conditions);
Vesicles_Released_AD = cell(1, num_conditions);

params_init_vr = [-0.04734 6.88 1.945 2.416];
lb_vr = [];
ub_vr = [];

params_vr_WT = cell(1, num_conditions);
params_vr_AD = cell(1, num_conditions);
Vesicles_Released_WT_fit = cell(1, num_conditions);
Vesicles_Released_AD_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end


    Vesicles_Released_WT{condition_index} = total_rel_ves_WT{condition_index}(100000, :);
    
   
    x_WT = log10(channel_number);
    Y_WT = Vesicles_Released_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)DoseResponseFit(params_WT, x_WT),params_init_vr, x_WT, Y_WT, lb_vr,ub_vr);
                                    
    params_vr_WT{condition_index} = params_WT;
    
    Vesicles_Released_WT_fit{condition_index} = DoseResponseFit(params_vr_WT{condition_index}, log10(channel_number));
    
  
end
   

%% Total Vesicle Released Synchronously (30ms after stimulation) vs Release Probability

sync_total_rel_ves_WT{condition_index} = cell(1, num_conditions);
sync_total_rel_ves_AD{condition_index} = cell(1, num_conditions);

SyncVesicles_Released_WT = cell(1, num_conditions);
SyncVesicles_Released_AD = cell(1, num_conditions);

params_init_vr = [-0.04734 6.88 1.945 2.416];
lb_vr = [];
ub_vr = [];

params_vr_WT = cell(1, num_conditions);
params_vr_AD = cell(1, num_conditions);
SyncVesicles_Released_WT_fit = cell(1, num_conditions);
SyncVesicles_Released_AD_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end

    
    sync_total_rel_ves_WT{condition_index} = cumtrapz(t, Sync_RelRate_WT{condition_index}, 1);
    SyncVesicles_Released_WT{condition_index} = sync_total_rel_ves_WT{condition_index}(end, :);    
   
    x_WT = log10(channel_number);
    Y_WT = SyncVesicles_Released_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)DoseResponseFit(params_WT, x_WT),params_init_vr, x_WT, Y_WT, lb_vr,ub_vr);
                                    
    params_vr_WT{condition_index} = params_WT;
    
    SyncVesicles_Released_WT_fit{condition_index} = DoseResponseFit(params_vr_WT{condition_index}, log10(channel_number));
    
  
end
   

%% Total Vesicle Released Asynchronously (30ms after stimulation) vs Release Probability

Async_total_rel_ves_WT{condition_index} = cell(1, num_conditions);
Async_total_rel_ves_AD{condition_index} = cell(1, num_conditions);

AsyncVesicles_Released_WT = cell(1, num_conditions);
AsyncVesicles_Released_AD = cell(1, num_conditions);

params_init_vr = [8.244 6.897 2.657 2.384 -2.134 -0.2979];
lb_vr = [];
ub_vr = [];

params_vr_WT = cell(1, num_conditions);
params_vr_AD = cell(1, num_conditions);
AsyncVesicles_Released_WT_fit = cell(1, num_conditions);
AsyncVesicles_Released_AD_fit = cell(1, num_conditions);

for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
       
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end

    
    Async_total_rel_ves_WT{condition_index} = cumtrapz(t, Async_RelRate_WT{condition_index}, 1);
    AsyncVesicles_Released_WT{condition_index} = Async_total_rel_ves_WT{condition_index}(end, :);
    
   
    x_WT = log10(channel_number);
    Y_WT = AsyncVesicles_Released_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)ExpDoseResponseFit(params_WT, x_WT),params_init_vr, x_WT, Y_WT, lb_vr,ub_vr);
                                    
    params_vr_WT{condition_index} = params_WT;
    
    AsyncVesicles_Released_WT_fit{condition_index} = ExpDoseResponseFit(params_vr_WT{condition_index}, log10(channel_number));
    
  
end


%% Total Vesicle Released Spontaneously (30ms after stimulation) vs Release Probability

Spont_total_rel_ves_WT{condition_index} = cell(1, num_conditions);
Spont_total_rel_ves_AD{condition_index} = cell(1, num_conditions);

SpontVesicles_Released_WT = cell(1, num_conditions);
SpontVesicles_Released_AD = cell(1, num_conditions);

params_init_vr = [0.009385 0.01337 55.8 0.0143];
lb_vr = [];
ub_vr = [];

params_vr_WT = cell(1, num_conditions);
params_vr_AD = cell(1, num_conditions);
SpontVesicles_Released_WT_fit = cell(1, num_conditions);
SpontVesicles_Released_AD_fit = cell(1, num_conditions);

for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
       
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end

    
    Spont_total_rel_ves_WT{condition_index} = cumtrapz(t, Spont_RelRate_WT{condition_index}, 1);
    SpontVesicles_Released_WT{condition_index} = Spont_total_rel_ves_WT{condition_index}(end, :);    
   
    x_WT = channel_number;
    Y_WT = SpontVesicles_Released_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)RevDoseResponseFit(params_WT, x_WT),params_init_vr, x_WT, Y_WT, lb_vr,ub_vr);
                                    
    params_vr_WT{condition_index} = params_WT;
    
    SpontVesicles_Released_WT_fit{condition_index} = RevDoseResponseFit(params_vr_WT{condition_index}, channel_number);
    
  
end


%% Time to Peak level (Time_To_Peak_Rate) Measured as time from simulation start to peak spike (release rate peak)


time_to_peak_WT = cell(1, num_conditions);
time_to_peak_AD = cell(1, num_conditions);

% set threshold release rate value as generic baseline rate
threshold = 0.01; % vesicles/ms
end_idx = length(t);


params_init_ttp = [197.3 1.075 2.931];
lb_ttp = [];
ub_ttp = [];

params_ttp_WT = cell(1, num_conditions);
params_ttp_AD = cell(1, num_conditions);
time_to_peak_WT_fit = cell(1, num_conditions);
time_to_peak_AD_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end
    
    [max_value_WT, max_index_WT] = max(rel_rate_WT{condition_index}, [], 1);

    for k=5:5:150
        
        j = int32(((k - 5)/5) + 1); 
        peak_index_WT = max_index_WT(1, j);
        time_to_peak_WT{condition_index}(1, j) = t(peak_index_WT);
      
            
    end
    
    %{
    x_WT = Pr_matrix_WT{condition_index};
    Y_WT = time_to_peak_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)ExponentialPolyFit(params_WT, x_WT),params_init_ttp, x_WT, Y_WT, lb_ttp,ub_ttp);
                                
    params_ttp_WT{condition_index} = params_WT;
    
    time_to_peak_WT_fit{condition_index} = ExponentialPolyFit(params_ttp_WT{condition_index}, Pr_matrix_WT{condition_index});
       
    %}
    
end

     
%% Cumulative AZ Calcium Concentration after stimulation Versus Channel Number 


params_init_cumm = [0.3937 17.53];
lb_cumm = [];
ub_cumm = [];

params_cumm_WT = cell(1, num_conditions);
params_cumm_AD = cell(1, num_conditions);
Cummul_VGCC_Ca_WT_fit = cell(1, num_conditions);
Cummul_VGCC_Ca_AD_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end


    x_WT = channel_number;
    Y_WT = Cummulative_VGCC_Ca_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)LinearFit(params_WT, x_WT),params_init_cumm, x_WT, Y_WT, lb_cumm,ub_cumm);
                                    
    params_cumm_WT{condition_index} = params_WT;
    
    Cummul_VGCC_Ca_WT_fit{condition_index} = LinearFit(params_cumm_WT{condition_index}, channel_number);
    
end

%% Time to basal level (Time_To_base) Measured as time from AP spike (release rate peak) to basal release


time_to_base_WT = cell(1, num_conditions);
time_to_base_AD = cell(1, num_conditions);

% set threshold release rate value as generic baseline rate
threshold = 0.001; % vesicles/ms
end_idx = length(t);


params_init_ttb = [197.3 1.075 2.931];
lb_ttb = [];
ub_ttb = [];

params_ttb_WT = cell(1, num_conditions);
params_ttb_AD = cell(1, num_conditions);
time_to_base_WT_fit = cell(1, num_conditions);
time_to_base_AD_fit = cell(1, num_conditions);


P2B_Cummulative_VGCC_Ca_WT = cell(1, num_conditions);
P2B_Cummulative_VGCC_Ca_AD = cell(1, num_conditions);
P2B_Cum_CouplingFlux_WT = cell(1, num_conditions);
P2B_Cum_CouplingFlux_AD = cell(1, num_conditions);

params_init_p2bcum = [47.41 0.8743 2.362];
lb_p2bcum = [];
ub_p2bcum = [];

params_p2bcum_WT = cell(1, num_conditions);
params_p2bcum_AD = cell(1, num_conditions);
P2B_Cummulative_VGCC_Ca_WT_fit = cell(1, num_conditions);
P2B_Cummulative_VGCC_Ca_AD_fit = cell(1, num_conditions);


for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
       
    sprintf(condition)
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end

    
        
    
    [max_value_WT, max_index_WT] = max(rel_rate_WT{condition_index}, [], 1);
  
    for k=5:5:150
        
        j = int32(((k - 5)/5) + 1);
        
        release_rate_WT = rel_rate_WT{condition_index}(max_index_WT(1, j):end, j);
        time_WT = t(1:end_idx-max_index_WT(1, j)+1);
        basal_index_WT = find(release_rate_WT <= threshold, 1);
        time_to_base_WT{condition_index}(1, j) = time_WT(basal_index_WT);
        
        ca_end_idx_WT = find(t >= str2double(num2str(time_to_base_WT{condition_index}(1, j) +...
                     t(max_index_WT(1, j)))), 1);
                 
        VGCC_Ca_WT = cumtrapz(t(max_index_WT(1, j): ca_end_idx_WT),...
                     Ca_VGCC_WT{condition_index}(max_index_WT(1, j):ca_end_idx_WT));
                              
        P2B_Cummulative_VGCC_Ca_WT{condition_index}(1, j) = VGCC_Ca_WT(end);
       
       
    end
    
    
   
    x_WT = Pr_matrix_WT{condition_index};
    Y_WT = time_to_base_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)ExponentialPolyFit(params_WT, x_WT),params_init_ttb, x_WT, Y_WT, lb_ttb,ub_ttb);
                                    
    params_ttb_WT{condition_index} = params_WT;
    
    time_to_base_WT_fit{condition_index} = ExponentialPolyFit(params_ttb_WT{condition_index}, Pr_matrix_WT{condition_index});
   
 
   
    x_WT = Pr_matrix_WT{condition_index};
    Y_WT = P2B_Cummulative_VGCC_Ca_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)ExponentialPolyFit(params_WT, x_WT),params_init_p2bcum, x_WT,...
                    Y_WT, lb_p2bcum, ub_p2bcum);
                                    
    params_p2bcum_WT{condition_index} = params_WT;
    
    P2B_Cummulative_VGCC_Ca_WT_fit{condition_index} = ExponentialPolyFit(params_p2bcum_WT{condition_index},...
                                                        Pr_matrix_WT{condition_index});
    
end
 
%% Plot overall, synchronous, asynchronous, and recruitment rates for synapse with 35 VGCCs

figure

subplot(4, 2, 1)
plot(t(1: 30000), rel_rate_WT{1}((1: 30000), 7),"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), rel_rate_WT{2}((1: 30000), 7),"k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), rel_rate_WT{3}((1: 30000), 7),"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(t(1: 30000), rel_rate_WT{4}((1: 30000), 7),'-', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%legend({'A\beta', 'IP_{3}R-NC', 'A\beta & IP_{3}R-NC', 'A\beta & IP_{3}R-HC'},'Location', 'northeast', 'FontSize',6)
legend({'IP_{3}R-NC', 'A\beta & IP_{3}R-NC', 'A\beta & IP_{3}R-HC'},'Location', 'northeast', 'FontSize',6)
ylabel('Release Rate (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
set(gca,'XTickLabelMode','auto') 
title('(A)', 'FontSize', 7);
hold off

subplot(4, 2, 2)
plot(channel_number, Pr_matrix_WT{1},"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Pr_matrix_WT{2} ,"k.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Pr_matrix_WT{3} ,"r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, Pr_matrix_WT{4} ,'.', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Pr_matrix_WT_fit{1} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Pr_matrix_WT_fit{2}, "k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Pr_matrix_WT_fit{3}, "r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, Pr_matrix_WT_fit{4}, '-', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
%legend boxoff
ylabel('Pr','FontSize',4,'FontWeight','bold','Color','k', 'Position', [-13, 0.5])
%ylim([0 1])
xlabel('Number of VGCCs','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);set(gca,'XTickLabelMode','auto')
title('(B)', 'FontSize', 7);
hold off


subplot(4, 2, 3)
plot(t(1: 30000), total_RRV_WT{1}((1: 30000), 7),"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), total_RRV_WT{2}((1: 30000), 7),"k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), total_RRV_WT{3}((1: 30000), 7), "r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(t(1: 30000), total_RRV_WT{4}((1: 30000), 7), '-', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
ylabel('RRP (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
set(gca,'XTickLabelMode','auto') 
title('(C)', 'FontSize', 7);

subplot(4, 2, 4)
plot(channel_number, Cummulative_VGCC_Ca_WT{1} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Cummul_VGCC_Ca_WT_fit{1} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Cummulative_VGCC_Ca_WT{2}, "k.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Cummul_VGCC_Ca_WT_fit{2} ,"k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Cummulative_VGCC_Ca_WT{3}, "r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Cummul_VGCC_Ca_WT_fit{3} ,"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, Cummulative_VGCC_Ca_WT{4}, '.', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, Cummul_VGCC_Ca_WT_fit{4} , '-', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
ylabel('Cumulative [Ca^{2+}]_{AZ} (\muM-ms)','FontSize',4,'FontWeight','bold','Color','k')
%ylim([0 90])
xlabel('Number of VGCCs','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);set(gca,'XTickLabelMode','auto') 
title('(D)', 'FontSize', 7);
hold off


subplot(4, 2, 5)
plot(t(1: 30000), Async_RelRate_WT{1}((1: 30000), 7)*1e03,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), Async_RelRate_WT{2}((1: 30000), 7)*1e03,"k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), Async_RelRate_WT{3}((1: 30000), 7)*1e03, "r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(t(1: 30000), Async_RelRate_WT{4}((1: 30000), 7)*1e03, '-', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
ylabel('Asynchronous Rate (10^{-3} ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
set(gca,'XTickLabelMode','auto') 
title('(E)', 'FontSize', 7);
hold off

subplot(4, 2, 6)
plot(channel_number, AsyncVesicles_Released_WT{1},"b.-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, AsyncVesicles_Released_WT_fit{1} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, AsyncVesicles_Released_WT{2}, "k.-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, AsyncVesicles_Released_WT_fit{2} ,"k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, AsyncVesicles_Released_WT{3}, "r.-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, AsyncVesicles_Released_WT_fit{3} ,"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, AsyncVesicles_Released_WT{4}, '.', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, AsyncVesicles_Released_WT_fit{4} , '-', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
ylabel('Total Asynchronous Release','FontSize',4,'FontWeight','bold','Color','k')
%ylim([0 5])
xlabel('Number of VGCCs','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);set(gca,'XTickLabelMode','auto') 
title('(F)', 'FontSize', 7);
hold off


%{
subplot(4, 2, 5)
plot(t(1: 30000), Spont_RelRate_WT{1}((1: 30000), 7)*1e03,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), Spont_RelRate_WT{2}((1: 30000), 7)*1e03,"k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), Spont_RelRate_WT{3}((1: 30000), 7)*1e03, "r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(t(1: 30000), Spont_RelRate_WT{4}((1: 30000), 7)*1e03, '-', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
ylabel('Spontaneous Rate (10^{-3} ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
set(gca,'XTickLabelMode','auto') 
title('(F)', 'FontSize', 7);
hold off


subplot(4, 2, 6)
plot(channel_number, SpontVesicles_Released_WT{1}*1e03,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, SpontVesicles_Released_WT_fit{1}*1e03 ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, SpontVesicles_Released_WT{2}*1e03, "k.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, SpontVesicles_Released_WT_fit{2}*1e03,"k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, SpontVesicles_Released_WT{3}*1e03, "r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, SpontVesicles_Released_WT_fit{3}*1e03 ,"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, SpontVesicles_Released_WT{4}*1e03, '.', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(channel_number, SpontVesicles_Released_WT_fit{4}*1e03 , '-', "color", "#A2142F", 'LineWidth',  0.85, 'MarkerSize', 8)
ylabel('Total Spontaneous Release','FontSize',4,'FontWeight','bold','Color','k')
%ylim([0 5])
xlabel('Number of VGCCs','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);set(gca,'XTickLabelMode','auto') 
title('(D)', 'FontSize', 7);
hold off
%}

subplot(4, 2, 7)
plot(Pr_matrix_WT{1}, time_to_base_WT{1} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{1}, time_to_base_WT_fit{1} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, time_to_base_WT{2} ,"k.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, time_to_base_WT_fit{2} ,"k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{3}, time_to_base_WT{3} ,"r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{3}, time_to_base_WT_fit{3} ,"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(Pr_matrix_WT{4}, time_to_base_WT{4} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
%hold on
%plot(Pr_matrix_WT{4}, time_to_base_WT_fit{4} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
%hold on
%ylim([0 25])
ylabel('Time-to-basal-rate (ms)','FontSize',4,'FontWeight','bold','Color','k')
xlabel('Pr','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);set(gca,'XTickLabelMode','auto') 
title('(G)', 'FontSize', 7);
hold off

subplot(4, 2, 8)
plot(Pr_matrix_WT{1}, P2B_Cummulative_VGCC_Ca_WT{1} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{1}, P2B_Cummulative_VGCC_Ca_WT_fit{1} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, P2B_Cummulative_VGCC_Ca_WT{2} ,"k.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, P2B_Cummulative_VGCC_Ca_WT_fit{2} ,"k-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{3}, P2B_Cummulative_VGCC_Ca_WT{3} ,"r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{3}, P2B_Cummulative_VGCC_Ca_WT_fit{3} ,"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
%plot(Pr_matrix_WT{4}, P2B_Cummulative_VGCC_Ca_WT{4} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
%hold on
%plot(Pr_matrix_WT{4}, P2B_Cummulative_VGCC_Ca_WT_fit{4} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
%hold on
ylabel({'Peak to base Cumulative'; '[Ca^{2+}]_{AZ} (\muM-ms)'},'FontSize',4,'FontWeight','bold','Color','k')
xlabel('Pr','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);set(gca,'XTickLabelMode','auto') 
title('(H)', 'FontSize', 7);
hold off


%Get Current Figure (GCF) & Set image size before saving image
width = 5.34*2.54;  % cm 
height = 5.98*2.54; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -djpeg -r1000 Figure_1
saveas(gcf, 'Figure_1', 'fig')
hold off


movefile Figure_1.jpg ../../results/Figure_1
movefile Figure_1.fig ../../results/Figure_1


%% %%%%%%%%%%%%%%%%%%%%%%%% FITTING   FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y_real = myode_for1(a, x)
        %%
        yD = a(1) - exp(-a(2)*(x + a(3)));
        y_real= yD;   
end 


function y_real = myode_for2(a, x)
        %%
        yD = a(1) + exp(-a(2)*(x - a(3)));
        y_real= yD;   
end 

function y_real = PolyExponentialFit(a, x)
        %%
        yD = (a(1) - x.^(a(2) * x.^a(3)))./x;
        y_real= yD;   
end 

function y_real = LogitFit(a, x)
        %% 
        yD = (a(1)).*(log(a(2).*(x.^a(3))./(1 - x)));
        y_real= yD;   
end 


function y_real = LinearFit(a, x)
        %%
        yD = a(1).*x + a(2);
        y_real= yD;   
end 


function y_real = DoseResponseFit(a, x)
        
        %% 
        yD = a(1) + (a(2) - a(1))./(1 + 10.^((a(3) - x).*a(4)));
        y_real= yD;   
end 

function y_real = RevDoseResponseFit(a, x)
        
        %% 
        yD = a(1) - (a(2) - a(1))./(1 + 10.^((a(3) - x).*a(4)));
        y_real= yD;   
end 

function y_real = SingleExpontialFit(a, x)
        
        %% 
        yD = a(1) + exp(a(2).*x + a(3));
        y_real= yD;   
end 

function y_real = GaussianFit(a, x)
        
        %% 
        yD = a(1)*exp(-((x - a(2))/a(3)).^2) + a(4)*exp(-((x -  a(5))/a(6)).^2);
        y_real= yD;   
end 

function y_real = BiexponentialFit(a, x)
        
        %%
               
        %yD = a(1)*exp(-x./a(2)) + a(3)*exp(-x./a(4)) + a(5)*exp(-x./a(6)) + a(7);
        
        yD = a(1).*exp(-x./a(2)) + a(3).*exp(-x./a(4));
        
        y_real= yD;    
end 


function y_real = ExponentialPolyFit(a, x)
        
        %%
           
        yD = a(1).*(x.^a(2)).*exp(-x.*a(3));
        
        y_real= yD;    
end

function y_real = PolynomialFit(a, x)
        
        %%
           
        yD = a(1).*(x.^2) + a(2).*x + a(3);
        
        y_real= yD;    
end

function y_real = PowerFit(a, x)
        
        %%
           
        yD = a(1).*(x.^a(2)) + a(3);
        
        y_real= yD;    
end

function y_real = ExpDoseResponseFit(a, x)

        yD = (a(2) - a(1))./(1 + 10.^((a(3) - x).*a(4))) + (exp(a(5).*(x.^a(6))));

        y_real = yD;
end

