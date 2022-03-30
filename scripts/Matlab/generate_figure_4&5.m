
clear all; close all; 


%% Data root directory

mkdir ../../results Figure_4
mkdir ../../results Figure_5

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
j = int32(((num_channels - 5)/5) + 1);    
figure

% Facilitation Computed Using Release Probability vs Stimulus Number


    subplot(4, 2, 1)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT_fit{2}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    legend({'WT', 'AD'},'Location', 'northeast', 'FontSize',2)
    str = {'Higher AD Coupling'};
    text(3, 0.7, str, 'FontSize', 5,'Color','k')
    ylabel('Facilitation (n^{th}/1^{st}) (Peak Rate)','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6)
    title('(A)', 'FontSize', 7);
    hold off
    

    subplot(4, 2, 2)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT_fit{2}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Facilitation (n^{th}/1^{st}) (P_{r})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6)
    title('(B)', 'FontSize', 7);
    hold off


    subplot(4, 2, 3)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT_fit{2}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Peak Rate (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto') 
    title('(C)', 'FontSize', 7);
    hold off
    

    subplot(4, 2, 4)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT_fit{2}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Pr','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(D)', 'FontSize', 7);
    hold off

    
% Plot RRP over time

    subplot(4, 2, 5)    
    plot(t(1: end), Fast_RRV_WT{2}((1: end), j) + Slow_RRV_WT{2}((1: end), j),"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Fast_RRV_AD{2}((1: end), j) + Slow_RRV_AD{2}((1: end), j),"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel('RRP','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(E)', 'FontSize', 7);
    hold off
    
   

% Plot [Ca^{2+}]_{AZ} (\muM) for synapse with 35 VGCCs   
    subplot(4, 2, 6) 
    plot(t(1: end), Ca_VGCC_WT{2}((1: end), j),"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Ca_VGCC_AD{2}((1: end), j),"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel('[Ca^{2+}]_{AZ} (\muM)','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
    xlim([0 450])
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6) 
    title('(F)', 'FontSize', 7);
    hold off


% Peak Asynchronous Rate

    subplot(4, 2, 7)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Async_Peak_WT{2}{num_channels}*1e03,"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Async_Peak_AD{2}{num_channels}*1e03,"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), 200), Async_Peak_WT_fit{2}{num_channels}*1e03,"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), 200), Async_Peak_AD_fit{2}{num_channels}*1e03, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Peak Asynchronous Rate (10^{-3} (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off'); a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto') 
    title('(G)', 'FontSize', 7);
    hold off
    

% Peak Synchronous Rate

    subplot(4, 2, 8)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT{2}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_AD{2}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT_fit{2}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_AD_fit{2}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Peak Synchronous Rate (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off'); a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto') 
    title('(H)', 'FontSize', 7);
    hold off
    

% Save Image

%Get Current Figure (GCF) & Set image size before saving image
width = 5.34*2.54;  % cm 
height = 5.98*2.54; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Figure_4
saveas(gcf, 'Figure_4', 'fig')
hold off


movefile Figure_4.png ../../results/Figure_4
movefile Figure_4.fig ../../results/Figure_4

%% Correlation between Spike Train and Release Events


[peaks, spike_times] = findpeaks(stimulus, t, 'MinPeakDistance', 20);

WT_corrcoef = cell(1, num_conditions);
AD_corrcoef = cell(1, num_conditions);
phase_WT = cell(1, num_conditions);
phase_AD = cell(1, num_conditions);
synchrony_WT = cell(1, num_conditions);
synchrony_AD = cell(1, num_conditions);
relative_synchrony_change = cell(1, num_conditions);
init_proba_WT = cell(1, num_conditions);
init_proba_AD = cell(1, num_conditions);

params_init_synchrony = [2.25 -1.75 -0.5027];
params_synchrony = cell(1, num_conditions);
relative_synchrony_change_fit = cell(1, num_conditions);
lb_sync = [];
ub_sync = [];

for condition_index = 1:3
    
    if condition_index == 1
        condition = "SameCoupling";
    elseif condition_index == 2
        condition = "AD_HigherCoupling";
    elseif condition_index == 3
        condition = "WT_HigherCoupling" ;
    end
    
    for num_channels=5:5:150
        j = int32(((num_channels - 5)/5) + 1);
        [peak_rel_WT, peak_rel_time_WT] = findpeaks(rel_rate_WT{condition_index}(:, j), t, 'MinPeakDistance', 20);
        [peak_rel_AD, peak_rel_time_AD] = findpeaks(rel_rate_AD{condition_index}(:, j), t, 'MinPeakDistance', 20);
       
        if length(peak_rel_time_WT) >= length(spike_times)
            stop_ind_WT = length(spike_times);
        elseif length(peak_rel_time_WT) <= length(spike_times)
            stop_ind_WT = length(peak_rel_time_WT);
        end 
        
        for m=1:stop_ind_WT-1
            phase_WT{condition_index}(j, m) = (spike_times(m) - peak_rel_time_WT(m))/...
                                              (peak_rel_time_WT(m+1)-peak_rel_time_WT(m));
        end 
            
        if length(peak_rel_time_AD) >= length(spike_times)
            stop_ind_AD = length(spike_times);
        elseif length(peak_rel_time_AD) <= length(spike_times)
            stop_ind_AD = length(peak_rel_time_AD);
        end 
        
        for m=1:stop_ind_AD-1
            phase_AD{condition_index}(j, m) = (spike_times(m) - peak_rel_time_AD(m))/...
                                              (peak_rel_time_AD(m+1)-peak_rel_time_AD(m));
        end  

        
        synchrony_WT{condition_index}(j) = real(mean(exp(0 + 2i*pi*phase_WT{condition_index}(j, :)),2));
        synchrony_AD{condition_index}(j) = real(mean(exp(0 + 2i*pi*phase_AD{condition_index}(j, :)),2));
        
        init_proba_WT{condition_index}(j) = Pr_WT{condition_index}{num_channels}(1);
        init_proba_AD{condition_index}(j) = Pr_AD{condition_index}{num_channels}(1);
        
        
        relative_synchrony_change{condition_index}(j) = (synchrony_AD{condition_index}(j) - synchrony_WT{condition_index}(j));
        
        

        x_synchrony = init_proba_WT{condition_index};
        Y_synchrony = relative_synchrony_change{condition_index};
        [params_sync] = lsqcurvefit(@(params_synchrony, x_synchrony)PolynomialFit(params_synchrony, x_synchrony),...
                             params_init_synchrony, x_synchrony, Y_synchrony, lb_sync, ub_sync);

        params_synchrony{condition_index} = params_sync;

        relative_synchrony_change_fit{condition_index} = PolynomialFit(params_synchrony{condition_index},...
                                                            init_proba_WT{condition_index});
                                                  
    end
        
end

%% Figure 5

figure

    subplot(2, 2, 1)
    ColorSet = jet(length(5:5:150));
    set(gca, 'ColorOrder', ColorSet);
    hold all;
    for num_channels=5:5:150
        j = int32(((num_channels - 5)/5) + 1);
        plot(phase_WT{2}(j, :), 'LineWidth', 0.7, 'MarkerSize', 8)
        hold on
    end
    num_channels = 35;
    j = int32(((num_channels - 5)/5) + 1);
    plot(phase_WT{2}(j, :), ".k", 'LineWidth', 0.7, 'MarkerSize', 6)
    ylabel({'Phase (\phi)) (WT)'},'FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel'); 
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(A)', 'FontSize', 7);
    hold off



    subplot(2, 2, 2)
    ColorSet = jet(length(5:5:150));
    set(gca, 'ColorOrder', ColorSet);
    hold all;
    for num_channels=5:5:150
        j = int32(((num_channels - 5)/5) + 1);
        plot(phase_AD{2}(j, :), 'LineWidth', 0.7, 'MarkerSize', 8)
        hold on
    end
    num_channels = 35;
    j = int32(((num_channels - 5)/5) + 1);
    plot(phase_AD{2}(j, :), ".k", 'LineWidth', 0.7, 'MarkerSize', 6)
    set(gca, 'ColorOrder', ColorSet);
    set(gca, 'Colormap', ColorSet);
    cb = colorbar('Ticks',[0.1, 0.9],...
             'TickLabels',["0.1 P_{r}",...
                           "0.9 P_{r}"], 'Location','eastoutside');
    cb.Position = [cb.Position(1)+0.1 cb.Position(2) cb.Position(3)-0.005 cb.Position(4)];
    caxis([0 1])
    ylabel({'Phase (\phi) (AD)'},'FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel'); 
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(B)', 'FontSize', 7);
    hold off


    subplot(2, 2, 3)
    plot(init_proba_WT{2}, synchrony_WT{2}, '.-b', 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(init_proba_AD{2}, synchrony_AD{2}, '.-r', 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel({'Synchrony'},'FontSize',4,'FontWeight','bold','Color','k')
    xlabel('P_{r}','FontSize',4,'FontWeight','bold','Color','k')
    legend({'WT', 'AD'},'Location', 'southwest', 'FontSize',3)
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel'); 
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title("(C)")
    hold off


    subplot(2, 2, 4)
    num_channels = (5:5:150);
    plot(init_proba_WT{2}, relative_synchrony_change{2}, '.', 'color', '#D95319', 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(init_proba_WT{2}, relative_synchrony_change_fit{2}, '-', 'color', '#D95319', 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel({'\Delta Synchrony'},'FontSize',4,'FontWeight','bold','Color','k')
    xlabel('P_{r}','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel'); 
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title("(D)")
    hold off

    %Get Current Figure (GCF) & Set image size before saving image
    width = 5.34*2.54;  % cm 
    height = 3*2.54; % cm
    set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

    %Set the resolution of 1000dpi and save the plot in TIFF format 
    print -dpng -r1000 Figure_5
    saveas(gcf, 'Figure_5', 'fig')
    hold off


movefile Figure_5.png ../../results/Figure_5
movefile Figure_5.fig ../../results/Figure_5

%% Plot zoomed in version of AZ calcium to show desynchronization

figure  

% Plot [Ca^{2+}]_{AZ} (\muM) for synapse with 35 VGCCs

    plot(t(1: end), Ca_VGCC_WT{2}((1: end), j),"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Ca_VGCC_AD{2}((1: end), j),"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel('[Ca^{2+}]_{AZ} (\muM)','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
    xlim([200 400])
    ylim([0 5])
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6) 
    title('(B)', 'FontSize', 5);
    hold off

%Get Current Figure (GCF) & Set image size before saving image
width = 8.2;  % cm 
height = 5.6; % cm
set(gcf, 'PaperPosition', [0, 0, width / 1.54, height / 1.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Figure_3
hold off


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