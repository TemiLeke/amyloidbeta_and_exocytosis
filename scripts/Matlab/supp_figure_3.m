
clear all; close all; 


%% Data root directory


mkdir ../../results supp_Figure_3

root = strcat(fileparts(fileparts(pwd)), "\data\Train_data\Preprocessed\");


%%   PPR, Probability, Release Rates and other plots

t = importdata(strcat(root, "time.txt"));
stimulus = importdata(strcat(root, "Stimulus_Train_20.txt"));


num_conditions = 3;                     % condition 1 corresponds to SamCoupling, 2 - AD Higher Coupling, 3 - WT Higher Coupling
num_channels = int32(((150 - 5)/5) + 1);
labels = cell(num_channels, 1);
channel_number = zeros(1, num_channels);

for k =5:5:150
    j = int32(((k - 5)/5) + 1);
    channel_number(1, j) = k;
end

total_rel_ves_WT = cell(1, num_conditions);
total_RRV_WT = cell(1, num_conditions);
Fast_RRV_WT = cell(1, num_conditions);
Slow_RRV_WT = cell(1, num_conditions); 
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
Fast_RRV_AD = cell(1, num_conditions);
Slow_RRV_AD = cell(1, num_conditions); 
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
    

    
    total_rel_ves_WT{k} = importdata(strcat(root,"WT_Facil_VesiclesReleased_", condition, ".csv"));
    total_RRV_WT{k} = importdata(strcat(root,"WT_Facil_RRV_", condition, ".csv"));
    Fast_RRV_WT{k} = importdata(strcat(root,"WT_Facil_Fast_RRV_", condition, ".csv"));
    Slow_RRV_WT{k} = importdata(strcat(root,"WT_Facil_Slow_RRV_", condition, ".csv"));
    rel_rate_WT{k} = importdata(strcat(root,"WT_Facil_ReleaseRate_", condition, ".csv"));
    Sync_RelRate_WT{k} = importdata(strcat(root,"WT_Facil_SyncRelRate_", condition, ".csv"));
    Async_RelRate_WT{k} =  importdata(strcat(root,"WT_Facil_AsyncRelRate_", condition, ".csv"));
    Spont_RelRate_WT{k} = importdata(strcat(root,"WT_Facil_SpontRelRate_", condition, ".csv"));
    Docked_Rate_WT{k} = importdata(strcat(root,"WT_Facil_DockedPoolRelRate_", condition, ".csv"));
    Reserved_Rate_WT{k} = importdata(strcat(root,"WT_Facil_ReservedPoolRelRate_", condition, ".csv"));
    Slow_Rate_WT{k} = importdata(strcat(root,"WT_Facil_SlowRelRate_", condition, ".csv"));
    Fast_Rate_WT{k} = importdata(strcat(root,"WT_Facil_FastRelRate_", condition, ".csv"));
    Ca_VGCC_WT{k} = importdata(strcat(root,"WT_Facil_Ca_VGCC_", condition, ".csv"));
    Ca_CYTO_WT{k} = importdata(strcat(root,"WT_Facil_Ca_CYTO_", condition, ".csv"));
    Ca_IP3R_WT{k} = importdata(strcat(root,"WT_Facil_Ca_IP3_", condition, ".csv")); 

    
    
    total_rel_ves_AD{k} = importdata(strcat(root,"AD_Facil_VesiclesReleased_", condition, ".csv"));
    total_RRV_AD{k} = importdata(strcat(root,"AD_Facil_RRV_", condition, ".csv"));
    Fast_RRV_AD{k} = importdata(strcat(root,"AD_Facil_Fast_RRV_", condition, ".csv"));
    Slow_RRV_AD{k} = importdata(strcat(root,"AD_Facil_Slow_RRV_", condition, ".csv"));
    rel_rate_AD{k} = importdata(strcat(root,"AD_Facil_ReleaseRate_", condition, ".csv"));
    Sync_RelRate_AD{k} = importdata(strcat(root,"AD_Facil_SyncRelRate_", condition, ".csv"));
    Async_RelRate_AD{k} = importdata(strcat(root,"AD_Facil_AsyncRelRate_", condition, ".csv"));
    Spont_RelRate_AD{k} = importdata(strcat(root,"AD_Facil_SpontRelRate_", condition, ".csv"));
    Docked_Rate_AD{k} = importdata(strcat(root,"AD_Facil_DockedPoolRelRate_", condition, ".csv"));
    Reserved_Rate_AD{k} = importdata(strcat(root,"AD_Facil_ReservedPoolRelRate_", condition, ".csv"));
    Slow_Rate_AD{k} = importdata(strcat(root,"AD_Facil_SlowRelRate_", condition, ".csv"));
    Fast_Rate_AD{k} = importdata(strcat(root,"AD_Facil_FastRelRate_", condition, ".csv"));
    Ca_VGCC_AD{k} = importdata(strcat(root,"AD_Facil_Ca_VGCC_", condition, ".csv"));
    Ca_CYTO_AD{k} = importdata(strcat(root,"AD_Facil_Ca_CYTO_", condition, ".csv"));
    Ca_IP3R_AD{k} = importdata(strcat(root,"AD_Facil_Ca_IP3_", condition, ".csv"));
    

end


%% Extract the Peaks of the stimulus train


[peaks, index_of_peaks] = findpeaks(stimulus);

peaks = peaks(peaks >= -20);
index_of_peaks = index_of_peaks(peaks >= -20);


%% Release Probability versus stimulus number For all Coupling Conditions


Pr_matrix_WT = cell(1, num_conditions);
Pr_matrix_AD = cell(1, num_conditions);
Pr_WT = cell(1, num_conditions);
Pr_AD = cell(1, num_conditions);

Peak_Rate_WT = cell(1, num_conditions);
Sync_Peak_WT = cell(1, num_conditions);
Async_Peak_WT = cell(1, num_conditions);


Peak_Rate_AD = cell(1, num_conditions);
Sync_Peak_AD = cell(1, num_conditions);
Async_Peak_AD = cell(1, num_conditions);


Facil_Pr_WT = cell(1, num_conditions);
Facil_Pr_AD = cell(1, num_conditions); 
Facil_PeakRate_WT = cell(1, num_conditions);
Facil_PeakRate_AD = cell(1, num_conditions); 
Sync_Facil_WT = cell(1, num_conditions);
Async_Facil_WT = cell(1, num_conditions);
Sync_Facil_AD = cell(1, num_conditions);
Async_Facil_AD = cell(1, num_conditions);


params_pr_WT = cell(1, num_conditions);
params_pr_AD = cell(1, num_conditions);

Peak_Rate_WT_fit = cell(1, num_conditions);
Peak_Rate_AD_fit = cell(1, num_conditions);
Pr_matrix_WT_fit = cell(1, num_conditions);
Pr_matrix_AD_fit = cell(1, num_conditions);
Pr_WT_fit = cell(1, num_conditions);
Pr_AD_fit = cell(1, num_conditions);
Facil_Pr_WT_fit = cell(1, num_conditions);
Facil_Pr_AD_fit = cell(1, num_conditions); 
Facil_PeakRate_WT_fit = cell(1, num_conditions);
Facil_PeakRate_AD_fit = cell(1, num_conditions);
Sync_Peak_WT_fit = cell(1, num_conditions);
Sync_Peak_AD_fit = cell(1, num_conditions);
Async_Peak_WT_fit = cell(1, num_conditions);
Async_Peak_AD_fit = cell(1, num_conditions);
Sync_Facil_WT_fit = cell(1, num_conditions);
Sync_Facil_AD_fit = cell(1, num_conditions);
Async_Facil_WT_fit = cell(1, num_conditions);
Async_Facil_AD_fit = cell(1, num_conditions);



for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
    
    for k=5:5:150
        j = int32(((k - 5)/5) + 1);
        index_of_base = 36900;
        duration = 21690;
        for i=1:length(peaks)
            
            if i==1
                Pr_matrix_WT{condition_index}{k}{i} = cumtrapz(t(1:index_of_base),...
                                                      Slow_Rate_WT{condition_index}(1:index_of_base, j)./5 + ...
                                                      Fast_Rate_WT{condition_index}(1:index_of_base, j)./5, 1);
                Pr_matrix_AD{condition_index}{k}{i} = cumtrapz(t(1:index_of_base),...
                                                      Slow_Rate_AD{condition_index}(1:index_of_base, j)./5 + ...
                                                      Fast_Rate_AD{condition_index}(1:index_of_base, j)./5, 1);
                                                  
                Peak_Rate_WT{condition_index}{k}(i) = max(rel_rate_WT{condition_index}(1:index_of_base, j));
                Sync_Peak_WT{condition_index}{k}(i) = max(Sync_RelRate_WT{condition_index}(1:index_of_base, j));
                Async_Peak_WT{condition_index}{k}(i) = max(Async_RelRate_WT{condition_index}(1:index_of_base, j));
                
                Peak_Rate_AD{condition_index}{k}(i) = max(rel_rate_AD{condition_index}(1:index_of_base, j));
                Sync_Peak_AD{condition_index}{k}(i) = max(Sync_RelRate_AD{condition_index}(1:index_of_base, j));
                Async_Peak_AD{condition_index}{k}(i) = max(Async_RelRate_AD{condition_index}(1:index_of_base, j));
                
            elseif i==length(peaks)
                Pr_matrix_WT{condition_index}{k}{i} = cumtrapz(t(index_of_base:end),...
                                                      Slow_Rate_WT{condition_index}(index_of_base:end, j)./5 + ...
                                                      Fast_Rate_WT{condition_index}(index_of_base:end, j)./5, 1);
                Pr_matrix_AD{condition_index}{k}{i} = cumtrapz(t(index_of_base:end),...
                                                      Slow_Rate_AD{condition_index}(index_of_base:end, j)./5 + ...
                                                      Fast_Rate_AD{condition_index}(index_of_base:end, j)./5, 1);
                                                  
                Peak_Rate_WT{condition_index}{k}(i) = max(rel_rate_WT{condition_index}(index_of_base:end, j));
                Sync_Peak_WT{condition_index}{k}(i) = max(Sync_RelRate_WT{condition_index}(index_of_base:end, j));
                Async_Peak_WT{condition_index}{k}(i) = max(Async_RelRate_WT{condition_index}(index_of_base:end, j));
                
                
                Peak_Rate_AD{condition_index}{k}(i) = max(rel_rate_AD{condition_index}(index_of_base:end, j)); 
                Sync_Peak_AD{condition_index}{k}(i) = max(Sync_RelRate_AD{condition_index}(index_of_base:end, j));
                Async_Peak_AD{condition_index}{k}(i) = max(Async_RelRate_AD{condition_index}(index_of_base:end, j));
                
            else
                
                Pr_matrix_WT{condition_index}{k}{i} = cumtrapz(t(index_of_base:index_of_base + duration),...
                                              Slow_Rate_WT{condition_index}(index_of_base:index_of_base + duration, j)./5 +...
                                              Fast_Rate_WT{condition_index}(index_of_base:index_of_base + duration, j)./5, 1);
                Pr_matrix_AD{condition_index}{k}{i} = cumtrapz(t(index_of_base:index_of_base + duration),...
                                              Slow_Rate_AD{condition_index}(index_of_base:index_of_base + duration, j)./5 +...
                                              Fast_Rate_AD{condition_index}(index_of_base:index_of_base + duration, j)./5, 1);
            
                Peak_Rate_WT{condition_index}{k}(i) = max(rel_rate_WT{condition_index}(index_of_base:...
                                                      index_of_base + duration, j));
                Sync_Peak_WT{condition_index}{k}(i) = max(Sync_RelRate_WT{condition_index}(index_of_base:...
                                                      index_of_base + duration, j));                               
                Async_Peak_WT{condition_index}{k}(i) = max(Async_RelRate_WT{condition_index}(index_of_base:...
                                                      index_of_base + duration, j));                                  
                                                  
                Peak_Rate_AD{condition_index}{k}(i) = max(rel_rate_AD{condition_index}(index_of_base:...
                                                      index_of_base + duration, j)); 
                Sync_Peak_AD{condition_index}{k}(i) = max(Sync_RelRate_AD{condition_index}(index_of_base:...
                                                      index_of_base + duration, j));                               
                Async_Peak_AD{condition_index}{k}(i) = max(Async_RelRate_AD{condition_index}(index_of_base:...
                                                      index_of_base + duration, j)); 
                                                  
            end
            
            index_of_base = index_of_base + duration;
            Pr_WT{condition_index}{k}(i) = Pr_matrix_WT{condition_index}{k}{i}(end);
            Facil_Pr_WT{condition_index}{k}(i) = Pr_WT{condition_index}{k}(i)/Pr_WT{condition_index}{k}(1);
            Facil_PeakRate_WT{condition_index}{k}(i) = Peak_Rate_WT{condition_index}{k}(i)/...
                                                        Peak_Rate_WT{condition_index}{k}(1);
            Sync_Facil_WT{condition_index}{k}(i) = Sync_Peak_WT{condition_index}{k}(i)/...
                                                    Sync_Peak_WT{condition_index}{k}(1);
            Async_Facil_WT{condition_index}{k}(i) = Async_Peak_WT{condition_index}{k}(i)/...
                                                    Async_Peak_WT{condition_index}{k}(1);                                               
                                                    
            Pr_AD{condition_index}{k}(i) = Pr_matrix_AD{condition_index}{k}{i}(end);
            Facil_Pr_AD{condition_index}{k}(i) = Pr_AD{condition_index}{k}(i)/Pr_AD{condition_index}{k}(1);
            Facil_PeakRate_AD{condition_index}{k}(i) = Peak_Rate_AD{condition_index}{k}(i)/...
                                                        Peak_Rate_AD{condition_index}{k}(1);
                                                    
            Sync_Facil_AD{condition_index}{k}(i) = Sync_Peak_AD{condition_index}{k}(i)/...
                                                    Sync_Peak_AD{condition_index}{k}(1);
            Async_Facil_AD{condition_index}{k}(i) = Async_Peak_AD{condition_index}{k}(i)/...
                                                    Async_Peak_AD{condition_index}{k}(1);                                          
        end
        
%########## Peak Rate Fit ############################################

        params_init_pr = [1.684 1/2.464 0.372 1/0.1968];
        lb_pr = [];
        ub_pr = [];
        x = linspace(1, length(index_of_peaks), length(index_of_peaks));
        Y_WT = Peak_Rate_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Peak_Rate_WT_fit{condition_index}{k}= BiexponentialFit(params_pr_WT{condition_index}, x);


        Y_AD = Peak_Rate_AD{condition_index}{k};
        [params_AD] = lsqcurvefit(@(params_AD, x)BiexponentialFit(params_AD, x),params_init_pr, x, Y_AD, lb_pr,ub_pr);

        params_pr_AD{condition_index} = params_AD;

        Peak_Rate_AD_fit{condition_index}{k}= BiexponentialFit(params_pr_AD{condition_index}, x);

        
%########## Synchronous Peak Rate Fit ############################################

        params_init_pr = [1.684 1/2.464 0.372 1/0.1968];
        lb_pr = [];
        ub_pr = [];
        x = linspace(1, length(index_of_peaks), length(index_of_peaks));
        Y_WT = Sync_Peak_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Sync_Peak_WT_fit{condition_index}{k}= BiexponentialFit(params_pr_WT{condition_index}, x);


        Y_AD = Sync_Peak_AD{condition_index}{k};
        [params_AD] = lsqcurvefit(@(params_AD, x)BiexponentialFit(params_AD, x),params_init_pr, x, Y_AD, lb_pr,ub_pr);

        params_pr_AD{condition_index} = params_AD;

        Sync_Peak_AD_fit{condition_index}{k}= BiexponentialFit(params_pr_AD{condition_index}, x);
        
%########## Asynchronous Peak Rate Fit ############################################
        

        params_init_pr = [0.01042 1/0.1439 0.01494 1/0.6828];
        lb_pr = [];
        ub_pr = [];
        x = linspace(1, length(index_of_peaks), length(index_of_peaks));
        Y_WT = Async_Peak_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Async_Peak_WT_fit{condition_index}{k}= BiexponentialFit(params_pr_WT{condition_index}, ...
                                                 linspace(1, length(index_of_peaks), 200));

        params_init_pr = [0.1888 -1.516 4.442];
        lb_pr = [];
        ub_pr = [];
        Y_AD = Async_Peak_AD{condition_index}{k};
        [params_AD] = lsqcurvefit(@(params_AD, x)ExponentialPolyFit2(params_AD, x),params_init_pr, x, Y_AD, lb_pr,ub_pr);

        params_pr_AD{condition_index} = params_AD;

        Async_Peak_AD_fit{condition_index}{k}= ExponentialPolyFit2(params_pr_AD{condition_index}, ...
                                                 linspace(1, length(index_of_peaks), 200));
     


%########## Peak Rate Facilitation Fit #################################

        params_init_pr = [3.752 1/2.464 0.8289 1/0.1968];
        lb_pr = [];
        ub_pr = [];
        Y_WT = Facil_PeakRate_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Facil_PeakRate_WT_fit{condition_index}{k}= BiexponentialFit(params_pr_WT{condition_index}, x);


        Y_AD = Facil_PeakRate_AD{condition_index}{k};
        [params_AD] = lsqcurvefit(@(params_AD, x)BiexponentialFit(params_AD, x),params_init_pr, x, Y_AD, lb_pr,ub_pr);

        params_pr_AD{condition_index} = params_AD;

        Facil_PeakRate_AD_fit{condition_index}{k}= BiexponentialFit(params_pr_AD{condition_index}, x);


        
%########## Synchronouus Facilitation Fit #################################

        params_init_pr = [3.752 1/2.464 0.8289 1/0.1968];
        lb_pr = [];
        ub_pr = [];
        Y_WT = Sync_Facil_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Sync_Facil_WT_fit{condition_index}{k}= BiexponentialFit(params_pr_WT{condition_index}, x);


        Y_AD = Sync_Facil_AD{condition_index}{k};
        [params_AD] = lsqcurvefit(@(params_AD, x)BiexponentialFit(params_AD, x),params_init_pr, x, Y_AD, lb_pr,ub_pr);

        params_pr_AD{condition_index} = params_AD;

        Sync_Facil_AD_fit{condition_index}{k}= BiexponentialFit(params_pr_AD{condition_index}, x);
        
        
        
%########## Asynchronouus Facilitation Fit #################################



        params_init_pr = [7.261 1/0.1444 -10.44 1/0.684];
        lb_pr = [];
        ub_pr = [];
        Y_WT = Async_Facil_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Async_Facil_WT_fit{condition_index}{k}= BiexponentialFit(params_pr_WT{condition_index}, ...
                                                 linspace(1, length(index_of_peaks), 200));


        params_init_pr = [70.62 -1.517 4.443];
        lb_pr = [];
        ub_pr = [];
        Y_AD = Async_Facil_AD{condition_index}{k};
        [params_AD] = lsqcurvefit(@(params_AD, x)ExponentialPolyFit2(params_AD, x),params_init_pr, x, Y_AD, lb_pr,ub_pr);

        params_pr_AD{condition_index} = params_AD;

        Async_Facil_AD_fit{condition_index}{k}= ExponentialPolyFit2(params_pr_AD{condition_index}, ...
                                                 linspace(1, length(index_of_peaks), 200));
                
        
       
        
        
%########## Release Probability Fit #####################################


        params_init_pr = [0.1659 -0.1499];
        lb_pr = [];
        ub_pr = [];
        x = linspace(1, length(index_of_peaks)-1, length(index_of_peaks)-1);
        Y_WT = Pr_WT{condition_index}{k}(1, 1:end-1);
        [params_WT] = lsqcurvefit(@(params_WT, x)ExponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Pr_WT_fit{condition_index}{k} = ExponentialFit(params_pr_WT{condition_index},...
                                        linspace(1, length(index_of_peaks), length(index_of_peaks)));

        params_init_pr = [0.04465 1/0.2243 0.1298 1/0.1279];
        Y_AD = Pr_AD{condition_index}{k}(1, 1:end-1);
        [params_AD] = lsqcurvefit(@(params_AD, x)BiexponentialFit(params_AD, x),params_init_pr, x, Y_AD, lb_pr,ub_pr);

        params_pr_AD{condition_index} = params_AD;

        Pr_AD_fit{condition_index}{k} = BiexponentialFit(params_pr_AD{condition_index},...
                                        linspace(1, length(index_of_peaks), length(index_of_peaks)));

 
        

%########## Release Probability Facilitation Fit #################################

        params_init_pr = [3.752 1/2.464 0.8289 1/0.1968];
        lb_pr = [];
        ub_pr = [];
        Y_WT = Facil_Pr_WT{condition_index}{k}(1, 1:end-1);
        [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Facil_Pr_WT_fit{condition_index}{k}= BiexponentialFit(params_pr_WT{condition_index},...
                                        linspace(1, length(index_of_peaks), length(index_of_peaks)));

        params_init_pr = [0.265 1/0.2311  0.9087 1/0.1296];
        Y_AD = Facil_Pr_AD{condition_index}{k}(1, 1:end-1);
        [params_AD] = lsqcurvefit(@(params_AD, x)BiexponentialFit(params_AD, x),params_init_pr, x, Y_AD, lb_pr,ub_pr);

        params_pr_AD{condition_index} = params_AD;

        Facil_Pr_AD_fit{condition_index}{k}= BiexponentialFit(params_pr_AD{condition_index},...
                                        linspace(1, length(index_of_peaks), length(index_of_peaks)));

        
        
    end

    
end

   
%% Plot

num_channels = 35;
figure

% Facilitation Computed Using Release Probability vs Stimulus Number


    subplot(3, 2, 1)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT_fit{2}{num_channels},"b--", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT_fit{3}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_AD_fit{3}{num_channels}, "r--", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on 
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT{3}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_AD{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    legend({'WT-NC', 'WT-HC', 'AD-HC','AD-NC'},'Location', 'northeast', 'FontSize',6)
    ylabel('Facilitation (n^{th}/1^{st}) (Peak Rate)','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6)
    title('(A)', 'FontSize', 7);
    hold off
    

    subplot(3, 2, 2)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT_fit{2}{num_channels},"b--", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT{3}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT_fit{3}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_AD{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_AD_fit{3}{num_channels}, "r--", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Facilitation (n^{th}/1^{st}) (P_{r})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6)
    title('(B)', 'FontSize', 7);
    hold off


    subplot(3, 2, 3)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT_fit{2}{num_channels},"b--", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT{3}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT_fit{3}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_AD{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_AD_fit{3}{num_channels}, "r--", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Peak Rate (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto') 
    title('(C)', 'FontSize', 7);
    hold off
    

    subplot(3, 2, 4)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT_fit{2}{num_channels},"b--", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT{3}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT_fit{3}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_AD{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_AD_fit{3}{num_channels}, "r--", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Pr','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(D)', 'FontSize', 7);
    hold off
 


% Peak Asynchronous Rate

    subplot(3, 2, 5)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Async_Peak_WT{2}{num_channels}*1e03,"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), 200), Async_Peak_WT_fit{2}{num_channels}*1e03,"b--", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Async_Peak_WT{3}{num_channels}*1e03,"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), 200), Async_Peak_WT_fit{3}{num_channels}*1e03,"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Async_Peak_AD{2}{num_channels}*1e03,"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), 200), Async_Peak_AD_fit{2}{num_channels}*1e03, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Async_Peak_AD{3}{num_channels}*1e03,"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), 200), Async_Peak_AD_fit{3}{num_channels}*1e03, "r--", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Peak Asynchronous Rate (10^{-3} (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off'); a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto') 
    title('(E)', 'FontSize', 7);
    hold off
    

% Peak Synchronous Rate

    subplot(3, 2, 6)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT_fit{2}{num_channels},"b--", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT{3}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT_fit{3}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_AD{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_AD_fit{3}{num_channels}, "r--", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Peak Synchronous Rate (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off'); a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto') 
    title('(F)', 'FontSize', 7);
    hold off
    

% Save Image

%Get Current Figure (GCF) & Set image size before saving image
width = 5.34*2.54;  % cm 
height = 5.98*2.54; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Figure_3
saveas(gcf, 'Figure_3', 'fig')
hold off


movefile Figure_3.png ../../results/supp_Figure_3
movefile Figure_3.fig ../../results/supp_Figure_3


%% %%%%%%%%%%%%%%%%%%%%%%%% FITTING   FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y_real = ExponentialFit(a, x)
        %%
        yD = a(1)*exp(a(2)*x);
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

function y_real = ExponentialPolyFit2(a, x)
        
        %%
           
        yD = a(1).*(x.^a(2)).*exp(-1./x.*(a(3)));
        
        y_real= yD;    
end


function y_real = PolynomialFit(a, x)
        
        %%
           
        yD = a(1).*(x.^2) + a(2).*x + a(3);
        
        y_real= yD;    
end

