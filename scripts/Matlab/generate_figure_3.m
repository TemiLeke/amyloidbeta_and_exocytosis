
clear all; close all; 


%% Data root directory

mkdir ../../results/Figure_3/

root = strcat(fileparts(fileparts(pwd)), "/data/APTrain_data/");


%%   PPR, Probability, Release Rates and other plots

base_dir = strcat(fileparts(fileparts(pwd)), "/data/");
t = importdata(strcat(base_dir, "aptrain_time.txt"));
stimulus = importdata(strcat(base_dir, "Stimulus_Train_20APs.txt"));



num_conditions = 3;                      % condition 1 corresponds to SamCoupling, 2 - AD Higher Coupling, 3 - WT Higher Coupling

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
    Fast_RRV_WT{k} = zeros(length(t), num_channels);
    Slow_RRV_WT{k} = zeros(length(t), num_channels);
    rel_rate_WT{k} = zeros(length(t), num_channels);
    Sync_RelRate_WT{k} = zeros(length(t), num_channels);
    Async_RelRate_WT{k} =  zeros(length(t), num_channels);
    Spont_RelRate_WT{k} = zeros(length(t), num_channels);
    Docked_Rate_WT{k} = zeros(length(t), num_channels);
    Reserved_Rate_WT{k} = zeros(length(t), num_channels);
    Slow_Rate_WT{k} = zeros(length(t), num_channels);
    Fast_Rate_WT{k} = zeros(length(t), num_channels);
    Ca_VGCC_WT{k} = zeros(length(t), num_channels);
    Ca_CYTO_WT{k} = zeros(length(t), num_channels);
    Ca_IP3R_WT{k} = zeros(length(t), num_channels);

    
    for channel_number=5:5:150

        j = int32(((channel_number - 5)/5) + 1);
    
        total_rel_ves_WT{k}(:, j) = importdata(strcat(root, condition, '/', "Train_RelVes_", num2str(channel_number), "_VGCC", ".txt"));
        total_RRV_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_RRV_", num2str(channel_number), "_VGCC", ".txt"));
        Fast_RRV_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_Fast_RRV_", num2str(channel_number), "_VGCC", ".txt"));
        Slow_RRV_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_Slow_RRV_", num2str(channel_number), "_VGCC", ".txt"));
        rel_rate_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_RelRate_", num2str(channel_number), "_VGCC", ".txt"));
        Sync_RelRate_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_SyncRelRate_", num2str(channel_number), ".txt"));
        Async_RelRate_WT{k}(:, j) =  importdata(strcat(root, condition, '/',"Train_AsyncRelRate_", num2str(channel_number), ".txt"));
        Spont_RelRate_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_SpontRelRate_", num2str(channel_number), ".txt"));
        Docked_Rate_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_Rate_DockedPool_", num2str(channel_number), ".txt"));
        Reserved_Rate_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_Rate_ReservePool_", num2str(channel_number), ".txt"));
        Slow_Rate_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_SlowRelRate_", num2str(channel_number), ".txt"));
        Fast_Rate_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_FastRelRate_", num2str(channel_number), ".txt"));
        Ca_VGCC_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_VGCC_Calcium_", num2str(channel_number), "_VGCC", ".txt"));
        Ca_CYTO_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_CYTO_Calcium_", num2str(channel_number), "_VGCC", ".txt"));
        Ca_IP3R_WT{k}(:, j) = importdata(strcat(root, condition, '/',"Train_IP3_Calcium_", num2str(channel_number), "_VGCC", ".txt")); 

    end 
    

end


%% Extract the Peaks and local minima_loc of the stimulus train


[peaks, index_of_peaks] = findpeaks(stimulus);

peaks = peaks(peaks >= -20);
index_of_peaks = index_of_peaks(peaks >= -20);


minima_loc = find(islocalmin(stimulus)==1);

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



for coupling_cond=1:length(coupling_conditions)

    condition = coupling_conditions(coupling_cond);
    
    if condition == "ip3r_nc"
        condition_index = 1;
    elseif condition == "ip3r_nc_and_abeta"
        condition_index = 2;
    elseif condition == "ip3r_hc_and_abeta"
        condition_index = 3;
    end
    
    for k=5:5:150
        j = int32(((k - 5)/5) + 1);
        index_of_base = 16800;
        duration = 21303;
        for i=1:length(peaks)
            
            if i==1
                Pr_matrix_WT{condition_index}{k}{i} = cumtrapz(t(1:index_of_base),...
                                                      Slow_Rate_WT{condition_index}(1:minima_loc(1), j)./5 + ...
                                                      Fast_Rate_WT{condition_index}(1:minima_loc(1), j)./5, 1);

                Peak_Rate_WT{condition_index}{k}(i) = max(rel_rate_WT{condition_index}(1:minima_loc(1), j));
                Sync_Peak_WT{condition_index}{k}(i) = max(Sync_RelRate_WT{condition_index}(1:minima_loc(1), j));
                Async_Peak_WT{condition_index}{k}(i) = max(Async_RelRate_WT{condition_index}(1:minima_loc(1), j));


%{
            elseif i==length(peaks)
                Pr_matrix_WT{condition_index}{k}{i} = cumtrapz(t(index_of_base:end),...
                                                      Slow_Rate_WT{condition_index}(index_of_base:end, j)./5 + ...
                                                      Fast_Rate_WT{condition_index}(index_of_base:end, j)./5, 1);
                                                  
                Peak_Rate_WT{condition_index}{k}(i) = max(rel_rate_WT{condition_index}(index_of_base:end, j));
                Sync_Peak_WT{condition_index}{k}(i) = max(Sync_RelRate_WT{condition_index}(index_of_base:end, j));
                Async_Peak_WT{condition_index}{k}(i) = max(Async_RelRate_WT{condition_index}(index_of_base:end, j));
                
                   
%}             
            else
                
                Pr_matrix_WT{condition_index}{k}{i} = cumtrapz(t(minima_loc(i-1):minima_loc(i)),...
                                              Slow_Rate_WT{condition_index}(minima_loc(i-1):minima_loc(i), j)./5 +...
                                              Fast_Rate_WT{condition_index}(minima_loc(i-1):minima_loc(i), j)./5, 1);
            
                Peak_Rate_WT{condition_index}{k}(i) = max(rel_rate_WT{condition_index}(minima_loc(i-1):minima_loc(i), j));
                Sync_Peak_WT{condition_index}{k}(i) = max(Sync_RelRate_WT{condition_index}(minima_loc(i-1):minima_loc(i), j));                               
                Async_Peak_WT{condition_index}{k}(i) = max(Async_RelRate_WT{condition_index}(minima_loc(i-1):minima_loc(i), j));                                  
                                                 
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
                                                                                            
        end
        
%########## Peak Rate Fit ############################################

        params_init_pr = [0.5616 -0.1377 0.1443];
        lb_pr = [];
        ub_pr = [];
        x = linspace(1, length(index_of_peaks), length(index_of_peaks));
        Y_WT = Peak_Rate_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)ExponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Peak_Rate_WT_fit{condition_index}{k} = ExponentialFit(params_pr_WT{condition_index}, x);
        
%########## Synchronous Peak Rate Fit ############################################

        params_init_pr = [0.5616 -0.1377 0.1443];
        lb_pr = [];
        ub_pr = [];
        x = linspace(1, length(index_of_peaks), length(index_of_peaks));
        Y_WT = Sync_Peak_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)ExponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Sync_Peak_WT_fit{condition_index}{k}= ExponentialFit(params_pr_WT{condition_index}, x);

        
%########## Asynchronous Peak Rate Fit ############################################
        

        params_init_pr = [0.00792 1/0.07272 -0.009445 1/0.3529];
        lb_pr = [];
        ub_pr = [];
        x = linspace(1, length(index_of_peaks), length(index_of_peaks));
        Y_WT = Async_Peak_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Async_Peak_WT_fit{condition_index}{k}= BiexponentialFit(params_pr_WT{condition_index}, ...
                                                 linspace(1, length(index_of_peaks), 200));


%########## Peak Rate Facilitation Fit #################################

        params_init_pr = [0.8669 -0.1377 0.2227];
        lb_pr = [];
        ub_pr = [];
        Y_WT = Facil_PeakRate_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)ExponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Facil_PeakRate_WT_fit{condition_index}{k}= ExponentialFit(params_pr_WT{condition_index}, x);

        
%########## Synchronouus Facilitation Fit #################################

        params_init_pr = [0.8669 -0.1377 0.2227];
        lb_pr = [];
        ub_pr = [];
        Y_WT = Sync_Facil_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)ExponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Sync_Facil_WT_fit{condition_index}{k}= ExponentialFit(params_pr_WT{condition_index}, x);

        
%########## Asynchronouus Facilitation Fit #################################

        params_init_pr = [9.734 1/0.07311 -11.64 1/0.3542];
        lb_pr = [];
        ub_pr = [];
        Y_WT = Async_Facil_WT{condition_index}{k};
        [params_WT] = lsqcurvefit(@(params_WT, x)BiexponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Async_Facil_WT_fit{condition_index}{k}= BiexponentialFit(params_pr_WT{condition_index}, ...
                                                 linspace(1, length(index_of_peaks), 200));
        
        
%########## Release Probability Fit #####################################


        params_init_pr = [0.1183 -0.1108 0.02964];
        lb_pr = [];
        ub_pr = [];
        x = linspace(1, length(index_of_peaks)-1, length(index_of_peaks)-1);
        Y_WT = Pr_WT{condition_index}{k}(1, 2:end);
        [params_WT] = lsqcurvefit(@(params_WT, x)ExponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Pr_WT_fit{condition_index}{k} = ExponentialFit(params_pr_WT{condition_index},...
                                        linspace(1, length(index_of_peaks), length(index_of_peaks)));
        

%########## Release Probability Facilitation Fit #################################

        params_init_pr = [2.17 -0.1108 0.4867];
        lb_pr = [];
        ub_pr = [];
        Y_WT = Facil_Pr_WT{condition_index}{k}(1, 2:end);
        [params_WT] = lsqcurvefit(@(params_WT, x)ExponentialFit(params_WT, x),params_init_pr, x, Y_WT, lb_pr,ub_pr);

        params_pr_WT{condition_index} = params_WT;

        Facil_Pr_WT_fit{condition_index}{k}= ExponentialFit(params_pr_WT{condition_index},...
                                        linspace(1, length(index_of_peaks), length(index_of_peaks)));


    end

    
end

   
%% Plot

num_channels = 35;
j = int32(((num_channels - 5)/5) + 1);    
figure

% Facilitation Computed Using Release Probability vs Stimulus Number


    subplot(4, 2, 1)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT{1}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT{2}{num_channels},"k.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on    
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT_fit{1}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT_fit{2}{num_channels}, "k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_PeakRate_WT_fit{3}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    legend({'IP_{3}R-NC', 'A\beta & IP_{3}R-NC', 'A\beta & IP_{3}R-HC'},'Location', 'northeast', 'FontSize',5)
    %str = {'Higher AD Coupling'};
    %text(3, 0.7, str, 'FontSize', 5,'Color','k')
    ylabel('Facilitation (n^{th}/1^{st}) (Peak Rate)','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6)
    title('(A)', 'FontSize', 7);
    hold off
    

    subplot(4, 2, 2)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT{1}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT{2}{num_channels},"k.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT_fit{1}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT_fit{2}{num_channels}, "k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Facil_Pr_WT_fit{3}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Facilitation (n^{th}/1^{st}) (P_{r})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6)
    title('(B)', 'FontSize', 7);
    hold off


    subplot(4, 2, 3)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT{1}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT{2}{num_channels},"k.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT_fit{1}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT_fit{2}{num_channels},"k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Peak_Rate_WT_fit{3}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Peak Rate (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto') 
    title('(C)', 'FontSize', 7);
    hold off
    
%{
    subplot(4, 2, 4)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT{1}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT{2}{num_channels},"k.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT_fit{1}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT_fit{2}{num_channels},"k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Pr_WT_fit{3}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Pr','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(D)', 'FontSize', 7);
    hold off

%}
    
% Plot RRP over time

    subplot(4, 2, 4)    
    plot(t(1: end), Fast_RRV_WT{1}((1: end), j) + Slow_RRV_WT{1}((1: end), j),"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Fast_RRV_WT{2}((1: end), j) + Slow_RRV_WT{2}((1: end), j),"k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Fast_RRV_WT{3}((1: end), j) + Slow_RRV_WT{3}((1: end), j),"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel('RRP','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto')
    title('(E)', 'FontSize', 7);
    hold off

% Plot Asynchronous Rel Rate over time

    subplot(4, 2, 5) 
    plot(t(1: end), Async_RelRate_WT{1}((1: end), j)*1e03,"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Async_RelRate_WT{2}((1: end), j)*1e03,"k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Async_RelRate_WT{3}((1: end), j)*1e03,"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel('Asynchronous Rate (10^{-3} ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
    xlim([0 450])
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6) 
    title('(F)', 'FontSize', 7);
    hold off
   

% Plot Asynchronous Rel Rate over time

    subplot(4, 2, 6) 
    plot(t(1: end), Sync_RelRate_WT{1}((1: end), j),"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Sync_RelRate_WT{2}((1: end), j),"k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Sync_RelRate_WT{3}((1: end), j),"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel('Synchronous Rate (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
    xlim([0 450])
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6) 
    title('(F)', 'FontSize', 7);
    hold off

%{    
% Plot [Ca^{2+}]_{AZ} (\muM) for synapse with 35 VGCCs   
    subplot(4, 2, 6) 
    plot(t(1: end), Ca_VGCC_WT{1}((1: end), j),"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Ca_VGCC_WT{2}((1: end), j),"k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(t(1: end), Ca_VGCC_WT{3}((1: end), j),"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel('[Ca^{2+}]_{AZ} (\muM)','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Time (ms)','FontSize',4,'FontWeight','bold','Color','k')
    xlim([0 450])
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6) 
    title('(F)', 'FontSize', 7);
    hold off
%}

% Peak Asynchronous Rate

    subplot(4, 2, 7)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Async_Peak_WT{1}{num_channels}*1e03,"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Async_Peak_WT{2}{num_channels}*1e03,"k.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Async_Peak_WT{3}{num_channels}*1e03,"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), 200), Async_Peak_WT_fit{1}{num_channels}*1e03,"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), 200), Async_Peak_WT_fit{2}{num_channels}*1e03, "k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), 200), Async_Peak_WT_fit{3}{num_channels}*1e03,"r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    ylabel('Peak_{Async} (10^{-3} (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
    xlabel('Stimulus number','FontSize',4,'FontWeight','bold','Color','k')
    set(gca, 'box', 'off'); a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6);
    set(gca,'XTickLabelMode','auto') 
    title('(G)', 'FontSize', 7);
    hold off
    

% Peak Synchronous Rate

    subplot(4, 2, 8)
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT{1}{num_channels},"b.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT{2}{num_channels},"k.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT{3}{num_channels},"r.", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT_fit{1}{num_channels},"b-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT_fit{2}{num_channels},"k-", 'LineWidth', 0.85, 'MarkerSize', 8)
    hold on
    plot(linspace(1, length(index_of_peaks), length(index_of_peaks)), Sync_Peak_WT_fit{3}{num_channels}, "r-", 'LineWidth', 0.85, 'MarkerSize', 8)
    ylabel('Peak_{Sync} (ms^{-1})','FontSize',4,'FontWeight','bold','Color','k')
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
print -djpeg -r1000 Figure_3
saveas(gcf, 'Figure_3', 'fig')
hold off


movefile Figure_3.jpg ../../results/Figure_3
movefile Figure_3.fig ../../results/Figure_3

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
        condition = "Train_HigherCoupling" ;
    end
    
    for num_channels=5:5:150
        j = int32(((num_channels - 5)/5) + 1);
        [peak_rel_WT, peak_rel_time_WT] = findpeaks(rel_rate_WT{condition_index}(:, j), t, 'MinPeakDistance', 20);
%        [peak_rel_AD, peak_rel_time_AD] = findpeaks(rel_rate_AD{condition_index}(:, j), t, 'MinPeakDistance', 20);
       
        if length(peak_rel_time_WT) >= length(spike_times)
            stop_ind_WT = length(spike_times);
        elseif length(peak_rel_time_WT) <= length(spike_times)
            stop_ind_WT = length(peak_rel_time_WT);
        end 
        
        for m=1:stop_ind_WT-1
            phase_WT{condition_index}(j, m) = (spike_times(m) - peak_rel_time_WT(m))/...
                                              (peak_rel_time_WT(m+1)-peak_rel_time_WT(m));
        end 
            

        
        synchrony_WT{condition_index}(j) = real(mean(exp(0 + 2i*pi*phase_WT{condition_index}(j, :)),2));
        
        init_proba_WT{condition_index}(j) = Pr_WT{condition_index}{num_channels}(1);
        
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



%% Figure 4


figure

    subplot(2, 2, 1)
    ColorSet = jet(length(5:5:150));
    set(gca, 'ColorOrder', ColorSet);
    hold all;
    for num_channels=5:5:150
        j = int32(((num_channels - 5)/5) + 1);
        if j ~= 27
            plot(phase_WT{2}(j, :), 'LineWidth', 0.7, 'MarkerSize', 8)
            hold on
        end
    end
    num_channels = 35;
    j = int32(((num_channels - 5)/5) + 1);
    plot(phase_WT{2}(j, :), ".k", 'LineWidth', 0.7, 'MarkerSize', 8)
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
    plot(phase_AD{2}(j, :), ".k", 'LineWidth', 0.7, 'MarkerSize', 8)
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
    print -djpeg -r1000 Figure_5
    saveas(gcf, 'Figure_5', 'fig')
    hold off


movefile Figure_5.jpg ../../results/Figure_5
movefile Figure_5.fig ../../results/Figure_5

%}
%% Plot zoomed in version of AZ calcium to show desynchronization

%{
    num_channels = 35;
    j = int32(((num_channels - 5)/5) + 1);    

    figure  
    % Plot [Ca^{2+}]_{AZ} (\muM) for synapse with 35 VGCCs

    plot(t(1: end), Ca_VGCC_WT{2}((1: end), j),"b-", 'LineWidth', 0.85, 'MarkerSize', 3)
    hold on
    plot(t(1: end), Ca_VGCC_AD{2}((1: end), j),"r-", 'LineWidth', 0.85, 'MarkerSize', 3)
    hold on
    xlim([200 400])
    ylim([0 5])
    set(gca, 'box', 'off');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',6) 
    hold off

    %Get Current Figure (GCF) & Set image size before saving image
    width = 3.2;  % cm 
    height = 1.6; % cm
    set(gcf, 'PaperPosition', [0, 0, width / 1.54, height / 1.54])

    %Set the resolution of 1000dpi and save the plot in TIFF format 
    print -djpeg -r1000 Figure_4_zoom
    saveas(gcf, 'Figure_4_zoom', 'fig')
    hold off

    movefile Figure_4_zoom.jpg ../../results/Figure_4
    movefile Figure_4_zoom.fig ../../results/Figure_4
%}

%% %%%%%%%%%%%%%%%%%%%%%%%% FITTING   FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y_real = ExponentialFit(a, x)
        %%
        yD = a(1)*exp(a(2)*x) + a(3);
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

