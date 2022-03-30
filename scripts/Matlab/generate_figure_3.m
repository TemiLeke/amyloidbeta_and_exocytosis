
clear all; close all; 


%% Data and figure  directory

mkdir ../../results Figure_3

root = strcat(fileparts(fileparts(pwd)), "\data\Paired_Pulse_data\Preprocessed\");

%%   PPR, Probability, Release Rates and other plots

t = importdata(strcat(root, "time.txt"));

num_conditions = 3;                     % condition 1 corresponds to SamCoupling, 2 - AD Higher Coupling, 3 - WT Higher Coupling
num_channels = int32(((150 - 5)/5) + 1);
labels = cell(num_channels, 1);
channel_number = zeros(1, num_channels);

for k =5:5:150
    j = int32(((k - 5)/5) + 1);
    channel_number(1, j) = k;
end

rel_ves_stim1_WT = cell(1, num_conditions);
rel_ves_stim2_WT = cell(1, num_conditions);
rel_rate_WT = cell(1, num_conditions);
PPR_WT = cell(1, num_conditions);
rel_proba_stim1_WT = cell(1, num_conditions);
rel_proba_stim2_WT = cell(1, num_conditions);
RRV_stim1_WT = cell(1, num_conditions);
RRV_stim2_WT = cell(1, num_conditions);
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

rel_ves_stim1_AD = cell(1, num_conditions);
rel_ves_stim2_AD = cell(1, num_conditions);
rel_rate_AD = cell(1, num_conditions);
PPR_AD = cell(1, num_conditions);
rel_proba_stim1_AD = cell(1, num_conditions);
rel_proba_stim2_AD = cell(1, num_conditions);
RRV_stim1_AD = cell(1, num_conditions);
RRV_stim2_AD = cell(1, num_conditions);
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
    
    
    rel_ves_stim1_WT{k} = importdata(strcat(root, "WT_PPR_VesiclesReleased_Stim1_", condition, ".csv"));
    rel_ves_stim2_WT{k} = importdata(strcat(root, "WT_PPR_VesiclesReleased_Stim2_", condition, ".csv"));
    rel_rate_WT{k} = importdata(strcat(root, "WT_PPR_ReleaseRate_", condition, ".csv"));
    PPR_WT{k} = importdata(strcat(root, "WT_PPR_Ratio_", condition, ".csv"));
    rel_proba_stim1_WT{k} = importdata(strcat(root, "WT_PPR_ReleaseProba_Stim1_", condition, ".csv"));
    rel_proba_stim2_WT{k} = importdata(strcat(root, "WT_PPR_ReleaseProba_Stim2_", condition, ".csv"));
    RRV_stim1_WT{k} = importdata(strcat(root, "WT_PPR_RRV_Stim1_", condition, ".csv"));
    RRV_stim2_WT{k} = importdata(strcat(root, "WT_PPR_RRV_Stim2_", condition, ".csv"));
    Sync_RelRate_WT{k} = importdata(strcat(root, "WT_PPR_SyncRelRate_", condition, ".csv"));
    Async_RelRate_WT{k} = importdata(strcat(root, "WT_PPR_AsyncRelRate_", condition, ".csv"));
    Spont_RelRate_WT{k} = importdata(strcat(root, "WT_PPR_SpontRelRate_", condition, ".csv"));
    Docked_Rate_WT{k} = importdata(strcat(root, "WT_PPR_DockedPoolRelRate_", condition, ".csv"));
    Reserved_Rate_WT{k} = importdata(strcat(root, "WT_PPR_ReservedPoolRelRate_", condition, ".csv"));
    Slow_Rate_WT{k} = importdata(strcat(root, "WT_PPR_SlowRelRate_", condition, ".csv"));
    Fast_Rate_WT{k} = importdata(strcat(root, "WT_PPR_FastRelRate_", condition, ".csv"));
    Ca_VGCC_WT{k} = importdata(strcat(root, "WT_PPR_Ca_VGCC_", condition, ".csv"));
    Ca_CYTO_WT{k} = importdata(strcat(root, "WT_PPR_Ca_CYTO_", condition, ".csv"));
    Ca_IP3R_WT{k} = importdata(strcat(root, "WT_PPR_Ca_IP3_", condition, ".csv")); 

    rel_ves_stim1_AD{k} = importdata(strcat(root, "AD_PPR_VesiclesReleased_Stim1_", condition, ".csv"));
    rel_ves_stim2_AD{k} = importdata(strcat(root, "AD_PPR_VesiclesReleased_Stim2_", condition, ".csv"));
    rel_rate_AD{k} = importdata(strcat(root, "AD_PPR_ReleaseRate_", condition, ".csv"));
    PPR_AD{k} = importdata(strcat(root, "AD_PPR_Ratio_", condition, ".csv"));
    rel_proba_stim1_AD{k} = importdata(strcat(root, "AD_PPR_ReleaseProba_Stim1_", condition, ".csv"));
    rel_proba_stim2_AD{k} = importdata(strcat(root, "AD_PPR_ReleaseProba_Stim2_", condition, ".csv"));
    RRV_stim1_AD{k} = importdata(strcat(root, "AD_PPR_RRV_Stim1_", condition, ".csv"));
    RRV_stim2_AD{k} = importdata(strcat(root, "AD_PPR_RRV_Stim2_", condition, ".csv"));
    Sync_RelRate_AD{k} = importdata(strcat(root, "AD_PPR_SyncRelRate_", condition, ".csv"));
    Async_RelRate_AD{k} = importdata(strcat(root, "AD_PPR_AsyncRelRate_", condition, ".csv"));
    Spont_RelRate_AD{k} = importdata(strcat(root, "AD_PPR_SpontRelRate_", condition, ".csv"));
    Docked_Rate_AD{k} = importdata(strcat(root, "AD_PPR_DockedPoolRelRate_", condition, ".csv"));
    Reserved_Rate_AD{k} = importdata(strcat(root, "AD_PPR_ReservedPoolRelRate_", condition, ".csv"));
    Slow_Rate_AD{k} = importdata(strcat(root, "AD_PPR_SlowRelRate_", condition, ".csv"));
    Fast_Rate_AD{k} = importdata(strcat(root, "AD_PPR_FastRelRate_", condition, ".csv"));
    Ca_VGCC_AD{k} = importdata(strcat(root, "AD_PPR_Ca_VGCC_", condition, ".csv"));
    Ca_CYTO_AD{k} = importdata(strcat(root, "AD_PPR_Ca_CYTO_", condition, ".csv"));
    Ca_IP3R_AD{k} = importdata(strcat(root, "AD_PPR_Ca_IP3_", condition, ".csv")); 

end


%% Paired-Pulse Ration (30 ms after simulation each pulse) versus Baseline Probability For all Coupling Conditions


Pr_matrix_WT_stim1 = cell(1, num_conditions);
Pr_matrix_WT_stim2 = cell(1, num_conditions);

Pr_matrix_AD_stim1 = cell(1, num_conditions);
Pr_matrix_AD_stim2 = cell(1, num_conditions);


params_init_ppr = [1.024 1.13 1.319];
lb_ppr = [];
ub_ppr = [];

params_ppr_WT = cell(1, num_conditions);
params_ppr_AD = cell(1, num_conditions);

PPR_WT_fit = cell(1, num_conditions);
PPR_AD_fit = cell(1, num_conditions);


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_Highe rCoupling" ;
    end
        
    Pr_matrix_WT_stim1{condition_index} = rel_proba_stim1_WT{condition_index}(30000, :);
    Pr_matrix_WT_stim2{condition_index} = rel_proba_stim2_WT{condition_index}(30000, :);
    
    Pr_matrix_AD_stim1{condition_index} = rel_proba_stim1_AD{condition_index}(30000, :);
    Pr_matrix_AD_stim2{condition_index} = rel_proba_stim2_AD{condition_index}(30000, :);    
    
    
    x = Pr_matrix_WT_stim1{condition_index};
    Y_WT = PPR_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)PolyExponentialFit(params_WT, x),params_init_ppr, x, Y_WT, lb_ppr,ub_ppr);
                                    
    params_ppr_WT{condition_index} = params_WT;
    
    PPR_WT_fit{condition_index} = PolyExponentialFit(params_ppr_WT{condition_index}, Pr_matrix_WT_stim1{condition_index});
    
   
    
    x = Pr_matrix_AD_stim1{condition_index};
    Y_AD = PPR_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x)PolyExponentialFit(params_AD, x),params_init_ppr, x, Y_AD, lb_ppr,ub_ppr);
                                    
    params_ppr_AD{condition_index} = params_AD;
    
    PPR_AD_fit{condition_index} = PolyExponentialFit(params_ppr_AD{condition_index}, Pr_matrix_AD_stim1{condition_index});
    
end



%% PR2 against PR1


params_init_pr2r1 = [-1.064 1.07 0.00194]; 
lb_pr2r1 = [];
ub_pr2r1 = [];

params_pr2r1_WT = cell(1, num_conditions);
params_pr2r1_AD = cell(1, num_conditions);


Pr_matrix_WT_stim2_fit = cell(1, num_conditions);
Pr_matrix_AD_stim2_fit = cell(1, num_conditions);


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_Highe rCoupling" ;
    end
    
    
    x = Pr_matrix_WT_stim1{condition_index};
    Y_WT = Pr_matrix_WT_stim2{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)PolynomialFit(params_WT, x),params_init_pr2r1, x, Y_WT, lb_pr2r1,ub_pr2r1);
                                    
    params_pr2r1_WT{condition_index} = params_WT;
    
    Pr_matrix_WT_stim2_fit{condition_index} = PolynomialFit(params_pr2r1_WT{condition_index}, Pr_matrix_WT_stim1{condition_index});
    
   
    
    x = Pr_matrix_AD_stim1{condition_index};
    Y_AD = Pr_matrix_AD_stim2{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x)PolynomialFit(params_AD, x),params_init_pr2r1, x, Y_AD, lb_pr2r1,ub_pr2r1);
                                    
    params_pr2r1_AD{condition_index} = params_AD;
    
    Pr_matrix_AD_stim2_fit{condition_index} = PolynomialFit(params_pr2r1_AD{condition_index}, Pr_matrix_AD_stim1{condition_index});
    
end


%% Cummulative Calcium Concentration after second pulse (30 ms after simulation each pulse) versus Baseline Probability For all Coupling Conditions


Cumulative_Ca_VGCC_WT_stim1 = cell(1, num_conditions);
Cumulative_Ca_VGCC_AD_stim1 = cell(1, num_conditions);

Cumulative_Ca_VGCC_WT_stim2 = cell(1, num_conditions);
Cumulative_Ca_VGCC_AD_stim2 = cell(1, num_conditions);


params_init_cum2 = [18.21 6.911 0.2096 ];
lb_cum2 = [];
ub_cum2 = [];

params_cum2_WT = cell(1, num_conditions);
params_cum2_AD = cell(1, num_conditions);

Cumulative_Ca_VGCC_WT_stim1_fit = cell(1, num_conditions);
Cumulative_Ca_VGCC_AD_stim1_fit = cell(1, num_conditions);

Cumulative_Ca_VGCC_WT_stim2_fit = cell(1, num_conditions);
Cumulative_Ca_VGCC_AD_stim2_fit = cell(1, num_conditions);


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_Highe rCoupling" ;
    end
        
    % Cumulative Ca in response to first pulse
    
    VGCC_Ca_WT = cumtrapz(t(44000:73999), Ca_VGCC_WT{condition_index}(44000:73999, :), 1);
    VGCC_Ca_AD = cumtrapz(t(44000:73999), Ca_VGCC_AD{condition_index}(44000:73999, :), 1);
     
    Cumulative_Ca_VGCC_WT_stim2{condition_index} = VGCC_Ca_WT(end,:);
    Cumulative_Ca_VGCC_AD_stim2{condition_index} = VGCC_Ca_AD(end,:);
    
   
    x = Pr_matrix_WT_stim1{condition_index};
    Y_WT = Cumulative_Ca_VGCC_WT_stim2{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)LogitFit(params_WT, x),params_init_cum2, x, Y_WT, lb_cum2,ub_cum2);
                                    
    params_cum2_WT{condition_index} = params_WT;
    
    Cumulative_Ca_VGCC_WT_stim2_fit{condition_index} = LogitFit(params_cum2_WT{condition_index}, Pr_matrix_WT_stim1{condition_index});
    
   
    
    x = Pr_matrix_AD_stim1{condition_index};
    Y_AD = Cumulative_Ca_VGCC_AD_stim2{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x)LogitFit(params_AD, x),params_init_cum2, x, Y_AD, lb_cum2,ub_cum2);
                                    
    params_cum2_AD{condition_index} = params_AD;
    
    Cumulative_Ca_VGCC_AD_stim2_fit{condition_index} = LogitFit(params_cum2_AD{condition_index}, Pr_matrix_AD_stim1{condition_index});
    
    
    % Cumulative Ca in response to first pulse 
    
    
    VGCC_Ca_WT = cumtrapz(t(1:30000), Ca_VGCC_WT{condition_index}(1:30000, :), 1);
    VGCC_Ca_AD = cumtrapz(t(1:30000), Ca_VGCC_AD{condition_index}(1:30000, :), 1);
     
    Cumulative_Ca_VGCC_WT_stim1{condition_index} = VGCC_Ca_WT(end,:);
    Cumulative_Ca_VGCC_AD_stim1{condition_index} = VGCC_Ca_AD(end,:);
    
   
    x = Pr_matrix_WT_stim1{condition_index};
    Y_WT = Cumulative_Ca_VGCC_WT_stim1{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x)LogitFit(params_WT, x),params_init_cum2, x, Y_WT, lb_cum2,ub_cum2);
                                    
    params_cum2_WT{condition_index} = params_WT;
    
    Cumulative_Ca_VGCC_WT_stim1_fit{condition_index} = LogitFit(params_cum2_WT{condition_index}, Pr_matrix_WT_stim1{condition_index});
    
   
    
    x = Pr_matrix_AD_stim1{condition_index};
    Y_AD = Cumulative_Ca_VGCC_AD_stim1{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x)LogitFit(params_AD, x),params_init_cum2, x, Y_AD, lb_cum2,ub_cum2);
                                    
    params_cum2_AD{condition_index} = params_AD;
    
    Cumulative_Ca_VGCC_AD_stim1_fit{condition_index} = LogitFit(params_cum2_AD{condition_index}, Pr_matrix_AD_stim1{condition_index});
    
end


%% STIM 2 -- Time to basal level (Time_To_base) Measured as time from AP spike (release rate peak) to basal release


time_to_base_stim2_WT = cell(1, num_conditions);
time_to_base_stim2_AD = cell(1, num_conditions);

% set threshold release rate value as generic baseline rate
threshold = 0.01; % vesicles/ms
stim2_end_idx = 73999;

params_init_ttb = [197.3 1.075 2.931];
lb_ttb = [];
ub_ttb = [];

params_ttb_WT = cell(1, num_conditions);
params_ttb_AD = cell(1, num_conditions);

time_to_base_stim2_WT_fit = cell(1, num_conditions);
time_to_base_stim2_AD_fit = cell(1, num_conditions);


P2B_Cummulative_VGCC_Ca_stim2_WT = cell(1, num_conditions);
P2B_Cummulative_VGCC_Ca_stim2_AD = cell(1, num_conditions);

params_init_p2bcum = [47.41 0.8743 2.362];
lb_p2bcum = [];
ub_p2bcum = [];

params_p2bcum_WT = cell(1, num_conditions);
params_p2bcum_AD = cell(1, num_conditions);
P2B_Cummulative_VGCC_Ca_stim2_WT_fit = cell(1, num_conditions);
P2B_Cummulative_VGCC_Ca_stim2_AD_fit = cell(1, num_conditions);


for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
        
    
    [max_value_WT, max_index_WT] = max(rel_rate_WT{condition_index}(44000:stim2_end_idx , :), [], 1);
    max_index_WT = max_index_WT + 44000;
    [max_value_AD, max_index_AD] = max(rel_rate_AD{condition_index}(44000:stim2_end_idx , :), [], 1);
    max_index_AD = max_index_AD + 44000;
  
    for k=5:5:150
        
        j = int32(((k - 5)/5) + 1);
        
        release_rate_WT = rel_rate_WT{condition_index}(max_index_WT(1, j):stim2_end_idx , j);
        time_WT = t(1:stim2_end_idx - max_index_WT(1, j)+1);
        basal_index_WT = find(release_rate_WT <= threshold, 1);
        time_to_base_stim2_WT{condition_index}(1, j) = time_WT(basal_index_WT);
        
        ca_end_idx_WT = find(t(1:stim2_end_idx) >= str2double(num2str(time_to_base_stim2_WT{condition_index}(1, j) +...
                     t(max_index_WT(1, j)))), 1);
                 
        VGCC_Ca_WT = cumtrapz(t(max_index_WT(1, j): ca_end_idx_WT),...
                     Ca_VGCC_WT{condition_index}(max_index_WT(1, j):ca_end_idx_WT));
                              
        P2B_Cummulative_VGCC_Ca_stim2_WT{condition_index}(1, j) = VGCC_Ca_WT(end);
        
        
        release_rate_AD = rel_rate_AD{condition_index}(max_index_AD(1, j):stim2_end_idx , j);
        time_AD = t(1:stim2_end_idx-max_index_AD(1, j)+1);
        basal_index_AD = find(release_rate_AD <= threshold, 1);
        time_to_base_stim2_AD{condition_index}(1, j) = time_AD(basal_index_AD);
        
        ca_end_idx_AD = find(t(1:stim2_end_idx) >= str2double(num2str(time_to_base_stim2_AD{condition_index}(1, j) +...
                     t(max_index_AD(1, j)))), 1);
                 
        VGCC_Ca_AD = cumtrapz(t(max_index_AD(1, j): ca_end_idx_AD),...
                     Ca_VGCC_AD{condition_index}(max_index_AD(1, j):ca_end_idx_AD));
                              
        P2B_Cummulative_VGCC_Ca_stim2_AD{condition_index}(1, j) = VGCC_Ca_AD(end);

    
    end
    
    
   
    x_WT = Pr_matrix_WT_stim1{condition_index};
    Y_WT = time_to_base_stim2_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)ExponentialPolyFit(params_WT, x_WT),params_init_ttb, x_WT, Y_WT, lb_ttb,ub_ttb);
                                    
    params_ttb_WT{condition_index} = params_WT;
    
    time_to_base_stim2_WT_fit{condition_index} = ExponentialPolyFit(params_ttb_WT{condition_index}, Pr_matrix_WT_stim1{condition_index});
    
    
    x_AD = Pr_matrix_AD_stim1{condition_index};
    Y_AD = time_to_base_stim2_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x_AD)ExponentialPolyFit(params_AD, x_AD),params_init_ttb, x_AD, Y_AD, lb_ttb,ub_ttb);
                                    
    params_ttb_AD{condition_index} = params_AD;
    
    time_to_base_stim2_AD_fit{condition_index} = ExponentialPolyFit(params_ttb_AD{condition_index}, Pr_matrix_AD_stim1{condition_index});
    
    
        
    
    x_WT = Pr_matrix_WT_stim1{condition_index};
    Y_WT = P2B_Cummulative_VGCC_Ca_stim2_WT{condition_index};
    [params_WT] = lsqcurvefit(@(params_WT, x_WT)ExponentialPolyFit(params_WT, x_WT),params_init_p2bcum, x_WT,...
                    Y_WT, lb_p2bcum, ub_p2bcum);
                                    
    params_p2bcum_WT{condition_index} = params_WT;
    
    P2B_Cummulative_VGCC_Ca_stim2_WT_fit{condition_index} = ExponentialPolyFit(params_p2bcum_WT{condition_index},...
                                                        Pr_matrix_WT_stim1{condition_index});
    
    
    x_AD = Pr_matrix_AD_stim1{condition_index};
    Y_AD = P2B_Cummulative_VGCC_Ca_stim2_AD{condition_index};
    [params_AD] = lsqcurvefit(@(params_AD, x_AD) ExponentialPolyFit(params_AD, x_AD),params_init_p2bcum, x_AD,...
                 Y_AD, lb_p2bcum, ub_p2bcum);
                                    
    params_p2bcum_AD{condition_index} = params_AD;
    
    P2B_Cummulative_VGCC_Ca_stim2_AD_fit{condition_index} =  ExponentialPolyFit(params_p2bcum_AD{condition_index},...
                                                      Pr_matrix_AD_stim1{condition_index});
    
    
    
end
                                   




%% Plot Images

num_channels = 35;
    
figure

    % Facilitation Computed Using Release Probability vs Stimulus Number
    subplot(3, 2, 1)
    plot(t(1: 80000), rel_rate_WT{2}((1: 80000), 7),"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: 80000), rel_rate_AD{2}((1: 80000), 7),"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    str = {'Higher AD Coupling'};
    text(6, 0.35, str, 'FontSize', 5, 'Color','k')
    ylabel('Release Rate (vesicles ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto') 
    title('(A)', 'FontSize', 7);
    hold off
    
 
    subplot(3, 2, 2)
    plot(Pr_matrix_WT_stim1{2}, PPR_WT{2} ,"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, PPR_AD{2} ,"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_WT_stim1{2}, PPR_WT_fit{2} ,"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, PPR_AD_fit{2}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    legend({'WT', 'AD'},'Location', 'northeast', 'FontSize',3)
    ylabel('PPR','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Pr','FontSize',4,'FontWeight','bold','Color','k')
    xlim([0 0.6])
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(B)', 'FontSize', 7);
    hold off


    subplot(3, 2, 3)
    plot(Pr_matrix_WT_stim1{2}, Cumulative_Ca_VGCC_WT_stim2{2} ,"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, Cumulative_Ca_VGCC_AD_stim2{2} ,"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_WT_stim1{2}, Cumulative_Ca_VGCC_WT_stim2_fit{2} ,"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, Cumulative_Ca_VGCC_AD_stim2_fit{2}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Cumulative [Ca^{2+}]_{AZ} (\muM-ms)','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Pr','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel'); 
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(C)', 'FontSize', 7);
    hold off    

    
    subplot(3, 2, 4)
    plot(Pr_matrix_WT_stim1{2}, Pr_matrix_WT_stim2{2} ,"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, Pr_matrix_AD_stim2{2} ,"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_WT_stim1{2}, Pr_matrix_WT_stim2_fit{2} ,"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, Pr_matrix_AD_stim2_fit{2}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Pr_{2}', 'FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Pr_{1}','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(D)', 'FontSize', 7);
    hold off

    
    subplot(3, 2, 5)
    plot(Pr_matrix_WT_stim1{2}, time_to_base_stim2_WT{2} ,"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_WT_stim1{2}, time_to_base_stim2_WT_fit{2} ,"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, time_to_base_stim2_AD{2}, "r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, time_to_base_stim2_AD_fit{2} ,"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Time-to-basal-rate_{2} (ms)','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Pr_{1}','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(E)', 'FontSize', 7);
    hold off



    subplot(3, 2, 6)
    plot(Pr_matrix_WT_stim1{2}, P2B_Cummulative_VGCC_Ca_stim2_WT{2} ,"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_WT_stim1{2}, P2B_Cummulative_VGCC_Ca_stim2_WT_fit{2} ,"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, P2B_Cummulative_VGCC_Ca_stim2_AD{2}, "r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(Pr_matrix_AD_stim1{2}, P2B_Cummulative_VGCC_Ca_stim2_AD_fit{2} ,"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel({'Peak to base Cumulative [Ca^{2+}]_{AZ}'; '(\muM-ms)'},'FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Pr_{1}','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel'); 
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(F)', 'FontSize', 7);
    hold off


    
%Get Current Figure (GCF) & Set image size before saving image
width = 5.34*2.54;  % cm 
height = 4.5*2.54; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Figure_3
saveas(gcf, 'Figure_3', 'fig')
hold off


movefile Figure_3.png ../../results/Figure_3
movefile Figure_3.fig ../../results/Figure_3


%% %%%%%%%%%%%%%%%%%%%%%%%% FITTING   FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y_real = myode_for1(a, x)
        
        %% Single Exponential Fit
        yD = a(1) - exp(-a(2)*(x + a(3)));
        y_real= yD;   
end 


function y_real = myode_for2(a, x)
        
        %% Single Exponential Fit
        yD = a(1) + exp(-a(2)*(x - a(3)));
        y_real= yD;   
end 

function y_real = PolyExponentialFit(a, x)
        
        %% Gaussian Fit
        yD = (a(1) - x.^(a(2) * x.^a(3)))./x;
        y_real= yD;   
end

function y_real = LogitFit(a, x)
        
        %% Gaussian Fit
        yD = (a(1)).*(log(a(2).*(x.^a(3))./(1 - x)));
        y_real= yD;   
end 


function y_real = LinearFit(a, x)
        
        %% Gaussian Fit
        yD = a(1).*x + a(2);
        y_real= yD;   
end 


function y_real = DoseResponseFit(a, x)
        
        %% Gaussian Fit
        yD = a(1) + (a(2) - a(1))./(1 + 10.^((a(3) - x).*a(4)));
        y_real= yD;   
end 


function y_real = SingleExpontialFit(a, x)
        
        %% Gaussian Fit
        yD = a(1) + exp(a(2).*x + a(3));
        y_real= yD;   
end 

function y_real = GaussianFit(a, x)
        
        %% Gaussian Fit
        yD = a(1)*exp(-((x - a(2))/a(3)).^2) + a(4)*exp(-((x -  a(5))/a(6)).^2);
        y_real= yD;   
end 

function y_real = BiexponentialFit(a, x)
        
        %% Double Exponential Fit
               
        %yD = a(1)*exp(-x./a(2)) + a(3)*exp(-x./a(4)) + a(5)*exp(-x./a(6)) + a(7);
        
        yD = a(1).*exp(-x./a(2)) + a(3).*exp(-x./a(4));
        
        y_real= yD;    
end 


function y_real = ExponentialPolyFit(a, x)
        
        %% Exponential with Polynomial Fit
           
        yD = a(1).*(x.^a(2)).*exp(-x.*a(3));
        
        y_real= yD;    
end

function y_real = PolynomialFit(a, x)
        
        %%
           
        yD = a(1).*(x.^2) + a(2).*x + a(3);
        
        y_real= yD;    
end
%}