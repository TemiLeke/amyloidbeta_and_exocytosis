clear all; close all; 


%% Data and figure  directory

mkdir ../../results supp_Figure_1

root = strcat(fileparts(fileparts(pwd)), "\data\SingleAP_data\Preprocessed\");

%%   PPR, Probability, Release Rates and other plots

t = importdata(strcat(root, "time.txt"));

num_conditions = 3;                      % condition 1 corresponds to SamCoupling, 2 - AD Higher Coupling, 3 - WT Higher Coupling
num_channels = int32(((150 - 5)/5) + 1);
labels = cell(num_channels, 1);
channel_number = zeros(1, num_channels);

for k =5:5:150
    j = int32(((k - 5)/5) + 1);
    channel_number(1, j) = k;
end

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



total_rel_ves_AD = cell(1, num_conditions);
total_RRV_AD = cell(1, num_conditions);
rel_proba_AD = cell(1, num_conditions);
rel_rate_AD = cell(1, num_conditions);
Sync_RelRate_AD = cell(1, num_conditions);
Async_RelRate_AD = cell(1, num_conditions);
Spont_RelRate_AD = cell(1, num_conditions);
Docked_Rate_AD = cell(1, num_conditions);
Reserved_Rate_AD = cell(1, num_conditions);
Slow_Rate_AD = cell(1, num_conditions);
Fast_Rate_AD = cell(1, num_conditions); 
Ca_VGCC_AD = cell(1, num_conditions);
Ca_CYTO_AD = cell(1, num_conditions);
Ca_IP3R_AD = cell(1, num_conditions);




for k = 1:num_conditions
    
    if k == 1
        condition = "SameCoupling";
    elseif k == 2
        condition = "AD_HigherCoupling";
    elseif k == 3
        condition = "WT_HigherCoupling" ;
    end
    

    
    total_rel_ves_WT{k} = importdata(strcat(root, "WT_VesiclesReleased_", condition, ".csv"));
    total_RRV_WT{k} = importdata(strcat(root, "WT_RRV_", condition, ".csv"));
    rel_proba_WT{k} = importdata(strcat(root, "WT_ReleaseProba_", condition, ".csv"));
    rel_rate_WT{k} = importdata(strcat(root, "WT_ReleaseRate_", condition, ".csv"));
    Sync_RelRate_WT{k} = importdata(strcat(root, "WT_SyncRelRate_", condition, ".csv"));
    Async_RelRate_WT{k} = importdata(strcat(root, "WT_AsyncRelRate_", condition, ".csv"));
    Spont_RelRate_WT{k} = importdata(strcat(root, "WT_SpontRelRate_", condition, ".csv"));
    Docked_Rate_WT{k} = importdata(strcat(root, "WT_DockedPoolRelRate_", condition, ".csv"));
    Reserved_Rate_WT{k} = importdata(strcat(root, "WT_ReservedPoolRelRate_", condition, ".csv"));
    Slow_Rate_WT{k} = importdata(strcat(root, "WT_SlowRelRate_", condition, ".csv"));
    Fast_Rate_WT{k} = importdata(strcat(root, "WT_FastRelRate_", condition, ".csv"));
    Ca_VGCC_WT{k} = importdata(strcat(root, "WT_Ca_VGCC_", condition, ".csv"));
    Ca_CYTO_WT{k} = importdata(strcat(root, "WT_Ca_CYTO_", condition, ".csv"));
    Ca_IP3R_WT{k} = importdata(strcat(root, "WT_Ca_IP3_", condition, ".csv")); 

    total_rel_ves_AD{k} = importdata(strcat(root, "AD_VesiclesReleased_", condition, ".csv"));
    total_RRV_AD{k} = importdata(strcat(root, "AD_RRV_", condition, ".csv"));
    rel_proba_AD{k} = importdata(strcat(root, "AD_ReleaseProba_", condition, ".csv"));
    rel_rate_AD{k} = importdata(strcat(root, "AD_ReleaseRate_", condition, ".csv"));
    Sync_RelRate_AD{k} = importdata(strcat(root, "AD_SyncRelRate_", condition, ".csv"));
    Async_RelRate_AD{k} = importdata(strcat(root, "AD_AsyncRelRate_", condition, ".csv"));
    Spont_RelRate_AD{k} = importdata(strcat(root, "AD_SpontRelRate_", condition, ".csv"));
    Docked_Rate_AD{k} = importdata(strcat(root, "AD_DockedPoolRelRate_", condition, ".csv"));
    Reserved_Rate_AD{k} = importdata(strcat(root, "AD_ReservedPoolRelRate_", condition, ".csv"));
    Slow_Rate_AD{k} = importdata(strcat(root, "AD_SlowRelRate_", condition, ".csv"));
    Fast_Rate_AD{k} = importdata(strcat(root, "AD_FastRelRate_", condition, ".csv"));
    Ca_VGCC_AD{k} = importdata(strcat(root, "AD_Ca_VGCC_", condition, ".csv"));
    Ca_CYTO_AD{k} = importdata(strcat(root, "AD_Ca_CYTO_", condition, ".csv"));
    Ca_IP3R_AD{k} = importdata(strcat(root, "AD_Ca_IP3_", condition, ".csv"));
    

end


%% Release Probability (30 ms after simulation) versus channel number For all Coupling Conditions

Pr_matrix_WT = cell(1, num_conditions);
Pr_matrix_AD = cell(1, num_conditions);

params_init_pr = [0.006876 1.061 1.862 2.326];
lb_pr = [];
ub_pr = [];

params_pr_WT = cell(1, num_conditions);
params_pr_AD = cell(1, num_conditions);
Pr_matrix_WT_fit = cell(1, num_conditions);
Pr_matrix_AD_fit = cell(1, num_conditions);


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
        
    Pr_matrix_WT{condition_index} = rel_proba_WT{condition_index}(30000, :);
    Pr_matrix_AD{condition_index} = rel_proba_AD{condition_index}(30000, :);

    
            
    %Pr_matrix_WT{condition_index} = total_rel_ves_WT{condition_index}(30000, :)./total_RRV_WT{condition_index}(30000-1, :);
    %Pr_matrix_AD{condition_index} = total_rel_ves_AD{condition_index}(30000, :)./total_RRV_AD{condition_index}(30000-1, :);
    
    x = log10(channel_number);
    Y_WT = Pr_matrix_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)DoseResponseFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);
                                    
    params_pr_WT{condition_index} = params_WT;
    
    Pr_matrix_WT_fit{condition_index} = DoseResponseFit(params_pr_WT{condition_index}, log10(channel_number));
    
    
    Y_AD = Pr_matrix_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x)DoseResponseFit(params_AD, x),params_init_pr, x, Y_AD, lb_pr,ub_pr);
                                    
    params_pr_AD{condition_index} = params_AD;
    
    Pr_matrix_AD_fit{condition_index} = DoseResponseFit(params_pr_AD{condition_index}, log10(channel_number));
    
end


%% Time to basal level (Time_To_base) Measured as time from AP spike (release rate peak) to basal release


time_to_base_WT = cell(1, num_conditions);
time_to_base_AD = cell(1, num_conditions);

% set threshold release rate value as generic baseline rate
threshold = 0.01; % vesicles/ms
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


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
        
    
    [max_value_WT, max_index_WT] = max(rel_rate_WT{condition_index}, [], 1);
    [max_value_AD, max_index_AD] = max(rel_rate_AD{condition_index}, [], 1);
  
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
       
        
        
        release_rate_AD = rel_rate_AD{condition_index}(max_index_AD(1, j):end, j);
        time_AD = t(1:end_idx-max_index_AD(1, j)+1);
        basal_index_AD = find(release_rate_AD <= threshold, 1);
        time_to_base_AD{condition_index}(1, j) = time_AD(basal_index_AD); 
        
        ca_end_idx_AD = find(t >= str2double(num2str(time_to_base_AD{condition_index}(1, j) +...
                     t(max_index_AD(1, j)))), 1);
                 
        VGCC_Ca_AD = cumtrapz(t(max_index_AD(1, j): ca_end_idx_AD),...
                     Ca_VGCC_AD{condition_index}(max_index_AD(1, j):ca_end_idx_AD));
        P2B_Cummulative_VGCC_Ca_AD{condition_index}(1, j) = VGCC_Ca_AD(end);

    
    end
    
    
   
    x_WT = Pr_matrix_WT{condition_index};
    Y_WT = time_to_base_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)ExponentialPolyFit(params_WT, x_WT),params_init_ttb, x_WT, Y_WT, lb_ttb,ub_ttb);
                                    
    params_ttb_WT{condition_index} = params_WT;
    
    time_to_base_WT_fit{condition_index} = ExponentialPolyFit(params_ttb_WT{condition_index}, Pr_matrix_WT{condition_index});
    
    
    x_AD = Pr_matrix_AD{condition_index};
    Y_AD = time_to_base_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x_AD)ExponentialPolyFit(params_AD, x_AD),params_init_ttb, x_AD, Y_AD, lb_ttb,ub_ttb);
                                    
    params_ttb_AD{condition_index} = params_AD;
    
    time_to_base_AD_fit{condition_index} = ExponentialPolyFit(params_ttb_AD{condition_index}, Pr_matrix_AD{condition_index});
    
    
 
    
    
    x_WT = Pr_matrix_WT{condition_index};
    Y_WT = P2B_Cummulative_VGCC_Ca_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)ExponentialPolyFit(params_WT, x_WT),params_init_p2bcum, x_WT,...
                    Y_WT, lb_p2bcum, ub_p2bcum);
                                    
    params_p2bcum_WT{condition_index} = params_WT;
    
    P2B_Cummulative_VGCC_Ca_WT_fit{condition_index} = ExponentialPolyFit(params_p2bcum_WT{condition_index},...
                                                        Pr_matrix_WT{condition_index});
    
    
    x_AD = Pr_matrix_AD{condition_index};
    Y_AD = P2B_Cummulative_VGCC_Ca_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x_AD) ExponentialPolyFit(params_AD, x_AD),params_init_p2bcum, x_AD,...
                 Y_AD, lb_p2bcum, ub_p2bcum);
                                    
    params_p2bcum_AD{condition_index} = params_AD;
    
    P2B_Cummulative_VGCC_Ca_AD_fit{condition_index} =  ExponentialPolyFit(params_p2bcum_AD{condition_index},...
                                                      Pr_matrix_AD{condition_index});
    
    
    
end
                      


%% Peak Release Rate (30 ms after simulation) versus channel number


Peak_RelRate_WT = cell(1, num_conditions);
Peak_RelRate_AD = cell(1, num_conditions);

params_init_prr = [-0.04734 6.88 1.945 2.416];
lb_prr = [];
ub_prr = [];

params_prr_WT = cell(1, num_conditions);
params_prr_AD = cell(1, num_conditions);
Peak_RelRate_WT_fit = cell(1, num_conditions);
Peak_RelRate_AD_fit = cell(1, num_conditions);


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
        
    Peak_RelRate_WT{condition_index} = max(rel_rate_WT{condition_index}, [], 1);
    Peak_RelRate_AD{condition_index} = max(rel_rate_AD{condition_index}, [], 1);
    
    
    x = log10(channel_number);
    Y_WT = Peak_RelRate_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)DoseResponseFit(params_WT, x),params_init_prr, x, Y_WT, lb_prr,ub_prr);
                                    
    params_prr_WT{condition_index} = params_WT;
    
    Peak_RelRate_WT_fit{condition_index} = DoseResponseFit(params_prr_WT{condition_index}, log10(channel_number));
    
    
    x = log10(channel_number);
    Y_AD = Peak_RelRate_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x)DoseResponseFit(params_AD, x),params_init_prr, x, Y_AD, lb_prr,ub_prr);
                                    
    params_prr_AD{condition_index} = params_AD;
    
    Peak_RelRate_AD_fit{condition_index} = DoseResponseFit(params_prr_AD{condition_index}, log10(channel_number));
    
end

%% Cumulative AZ Calcium Concentration (30 ms) after stimulation Versus Release Probability 

Cummulative_VGCC_Ca_WT = cell(1, num_conditions);
Cummulative_VGCC_Ca_AD = cell(1, num_conditions);


params_init_cum = [18.21 6.911 0.2096 ];
lb_cum = [];
ub_cum = [];

params_cum_WT = cell(1, num_conditions);
params_cum_AD = cell(1, num_conditions);
Cummulative_VGCC_Ca_WT_fit = cell(1, num_conditions);
Cummulative_VGCC_Ca_AD_fit = cell(1, num_conditions);


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
    
    VGCC_Ca_WT = cumtrapz(t, Ca_VGCC_WT{condition_index}, 1);
    VGCC_Ca_AD = cumtrapz(t, Ca_VGCC_AD{condition_index}, 1);
     
    Cummulative_VGCC_Ca_WT{condition_index} = VGCC_Ca_WT(30000,:);
    Cummulative_VGCC_Ca_AD{condition_index} = VGCC_Ca_AD(30000,:);
    
    x_WT = Pr_matrix_WT{condition_index};
    Y_WT = Cummulative_VGCC_Ca_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)LogitFit(params_WT, x_WT),params_init_cum, x_WT, Y_WT, lb_cum,ub_cum);
                                    
    params_cum_WT{condition_index} = params_WT;
    
    Cummulative_VGCC_Ca_WT_fit{condition_index} = LogitFit(params_cum_WT{condition_index},Pr_matrix_WT{condition_index});
    
    
    x_AD = Pr_matrix_AD{condition_index};
    Y_AD = Cummulative_VGCC_Ca_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x_AD)LogitFit(params_AD, x_AD),params_init_cum, x_AD, Y_AD, lb_cum,ub_cum);
                                    
    params_cum_AD{condition_index} = params_AD;
    
    Cummulative_VGCC_Ca_AD_fit{condition_index} = LogitFit(params_cum_AD{condition_index}, Pr_matrix_AD{condition_index});
    
end


%% Total Vesicle Released (30ms after stimulation) vs Release Probability


Vesicles_Released_WT = cell(1, num_conditions);
Vesicles_Released_AD = cell(1, num_conditions);

params_init_vr = [-0.04734 6.88 1.945 2.416];
lb_vr = [];
ub_vr = [];

params_vr_WT = cell(1, num_conditions);
params_vr_AD = cell(1, num_conditions);
Vesicles_Released_WT_fit = cell(1, num_conditions);
Vesicles_Released_AD_fit = cell(1, num_conditions);


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
        
    Vesicles_Released_WT{condition_index} = total_rel_ves_WT{condition_index}(30000, :);
    Vesicles_Released_AD{condition_index} = total_rel_ves_AD{condition_index}(30000, :);
    
   
    x_WT = log10(channel_number);
    Y_WT = Vesicles_Released_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)DoseResponseFit(params_WT, x_WT),params_init_vr, x_WT, Y_WT, lb_vr,ub_vr);
                                    
    params_vr_WT{condition_index} = params_WT;
    
    Vesicles_Released_WT_fit{condition_index} = DoseResponseFit(params_vr_WT{condition_index}, log10(channel_number));
    
    
    x_AD = log10(channel_number);
    Y_AD = Vesicles_Released_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x_AD)DoseResponseFit(params_AD, x_AD),params_init_vr, x_AD, Y_AD, lb_vr,ub_vr);
                                    
    params_vr_AD{condition_index} = params_AD;
    
    Vesicles_Released_AD_fit{condition_index} = DoseResponseFit(params_vr_AD{condition_index}, log10(channel_number));
    
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



for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
        
    
    [max_value_WT, max_index_WT] = max(rel_rate_WT{condition_index}, [], 1);
    [max_value_AD, max_index_AD] = max(rel_rate_AD{condition_index}, [], 1);
  
    for k=5:5:150
        
        j = int32(((k - 5)/5) + 1); 
        peak_index_WT = max_index_WT(1, j);
        time_to_peak_WT{condition_index}(1, j) = t(peak_index_WT);
        
        
        peak_index_AD = max_index_AD(1, j);
        time_to_peak_AD{condition_index}(1, j) = t(peak_index_AD);
            
    end
    
    %{
    x_WT = Pr_matrix_WT{condition_index};
    Y_WT = time_to_peak_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)ExponentialPolyFit(params_WT, x_WT),params_init_ttp, x_WT, Y_WT, lb_ttp,ub_ttp);
                                
    params_ttp_WT{condition_index} = params_WT;
    
    time_to_peak_WT_fit{condition_index} = ExponentialPolyFit(params_ttp_WT{condition_index}, Pr_matrix_WT{condition_index});
    
        
    x_AD = Pr_matrix_AD{condition_index};
    Y_AD = time_to_peak_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x_AD)ExponentialPolyFit(params_AD, x_AD),params_init_ttp, x_AD, Y_AD, lb_ttp,ub_ttp);
                                
    params_ttp_AD{condition_index} = params_AD;
    
    time_to_peak_AD_fit{condition_index} = ExponentialPolyFit(params_ttp_AD{condition_index}, Pr_matrix_AD{condition_index});
    
    %}
    
end
     

%% Cumulative AZ Calcium Concentration (30 ms) after stimulation Versus Channel Number 


params_init_cumm = [0.3937 17.53];
lb_cumm = [];
ub_cumm = [];

params_cumm_WT = cell(1, num_conditions);
params_cumm_AD = cell(1, num_conditions);
Cummul_VGCC_Ca_WT_fit = cell(1, num_conditions);
Cummul_VGCC_Ca_AD_fit = cell(1, num_conditions);


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
    
    
    x_WT = channel_number;
    Y_WT = Cummulative_VGCC_Ca_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)LinearFit(params_WT, x_WT),params_init_cumm, x_WT, Y_WT, lb_cumm,ub_cumm);
                                    
    params_cumm_WT{condition_index} = params_WT;
    
    Cummul_VGCC_Ca_WT_fit{condition_index} = LinearFit(params_cumm_WT{condition_index}, channel_number);
    
    
    x_AD = channel_number;
    Y_AD = Cummulative_VGCC_Ca_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x_AD)LinearFit(params_AD, x_AD),params_init_cumm, x_AD, Y_AD, lb_cumm,ub_cumm);
                                    
    params_cumm_AD{condition_index} = params_AD;
    
    Cummul_VGCC_Ca_AD_fit{condition_index} = LinearFit(params_cumm_AD{condition_index}, channel_number);
    
end
                                    
                                   

%% Plot release rate, release probability, TTB and Cummulative Calcium profile for synapse with 35 VGCCs

figure

subplot(2, 2, 1)
plot(t(1: 30000), rel_rate_WT{2}((1: 30000), 7),"b--", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), rel_rate_WT{3}((1: 30000), 7),"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), rel_rate_AD{2}((1: 30000), 7),"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(t(1: 30000), rel_rate_AD{3}((1: 30000), 7),"r--", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
legend({'WT-NC', 'WT-HC', 'AD-HC','AD-NC'},'Location', 'northeast', 'FontSize',6)
%str = {'Higher AD Coupling with Parameters;','AD: k = 15 , K_{c} = 10 (\muM)',...
%      'WT: k = 5 , K_{c} = 20 (\muM)'};
%text(6, 0.35, str, 'FontSize', 5, 'Color','k')
ylabel('Release Rate (vesicles ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
set(gca,'XTickLabelMode','auto') 
title('(A)', 'FontSize', 7);
hold off


subplot(2, 2, 2)
plot(channel_number, Vesicles_Released_WT{2} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Vesicles_Released_WT_fit{2} ,"b--", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Vesicles_Released_WT{3} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Vesicles_Released_WT_fit{3} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Vesicles_Released_AD{2}, "r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Vesicles_Released_AD_fit{2} ,"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Vesicles_Released_AD{3}, "r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(channel_number, Vesicles_Released_AD_fit{3} ,"r--", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
ylabel('Vesicles Released','FontSize',4,'FontWeight','bold','Color','k')
ylim([0 5])
xlabel('Number of VGCCs','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);set(gca,'XTickLabelMode','auto') 
title('(B)', 'FontSize', 7);
hold off


subplot(2, 2, 3)
plot(Pr_matrix_WT{2}, time_to_base_WT{2} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, time_to_base_WT_fit{2} ,"b--", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, time_to_base_WT{3} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, time_to_base_WT_fit{3} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_AD{2}, time_to_base_AD{2}, "r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_AD{2}, time_to_base_AD_fit{2} ,"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_AD{2}, time_to_base_AD{3}, "r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_AD{2}, time_to_base_AD_fit{3} ,"r--", 'LineWidth',  0.85, 'MarkerSize', 8)
ylabel('Time-to-basal-rate (ms)','FontSize',4,'FontWeight','bold','Color','k')
ylim([0 25])
xlabel('Pr','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);set(gca,'XTickLabelMode','auto') 
title('(C)', 'FontSize', 7);
hold off




subplot(2, 2, 4)
plot(Pr_matrix_WT{2}, P2B_Cummulative_VGCC_Ca_WT{2} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, P2B_Cummulative_VGCC_Ca_WT_fit{2} ,"b--", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, P2B_Cummulative_VGCC_Ca_WT{3} ,"b.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_WT{2}, P2B_Cummulative_VGCC_Ca_WT_fit{3} ,"b-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_AD{2}, P2B_Cummulative_VGCC_Ca_AD{2}, "r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_AD{2}, P2B_Cummulative_VGCC_Ca_AD_fit{2} ,"r-", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_AD{2}, P2B_Cummulative_VGCC_Ca_AD{3}, "r.", 'LineWidth',  0.85, 'MarkerSize', 8)
hold on
plot(Pr_matrix_AD{2}, P2B_Cummulative_VGCC_Ca_AD_fit{3} ,"r--", 'LineWidth',  0.85, 'MarkerSize', 8)
ylabel({'Peak to base Cumulative'; '[Ca^{2+}]_{AZ} (\muM-ms)'},'FontSize',4,'FontWeight','bold','Color','k')
xlabel('Pr','FontSize',4,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);set(gca,'XTickLabelMode','auto') 
title('(D)', 'FontSize', 7);
hold off

%Get Current Figure (GCF) & Set image size before saving image
width = 5.34*2.54;  % cm 
height = 5.98*2.54; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Figure_1
saveas(gcf, 'Figure_1', 'fig')
hold off


movefile Figure_1.png ../../results/supp_Figure_1
movefile Figure_1.fig ../../results/supp_Figure_1


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



%}