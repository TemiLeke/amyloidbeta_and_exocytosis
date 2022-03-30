
clear all; close all; 

%%  Load Data
mkdir ../../results Figure_1

load ../../data/Figure1_data/DoseResp_of_Calyx.csv      % Data obtained from Sun, J et al. (2007). https://doi.org/10.1038/nature06308
load ../../data/Figure1_data/TimeToPeak_of_Calyx.csv    % Data obtained from Sun, J et al. (2007). https://doi.org/10.1038/nature06308
load ../../data/Figure1_data/SpontaneousRelRate_from_Suhita_2010.csv % Data obtained from Nadkarni, S. et al. (2010) doi:10.1371/journal.pcbi.1000983
load ../../data/Figure1_data/DoseResp_of_cMFBs.csv      % Data obtained from Eshra et al. eLife (2021) DOI: https://doi.org/10.7554/eLife.704081
load ../../data/Figure1_data/TimeToPeak_of_cMFBs.csv    % Data obtained from Eshra et al. (2021) DOI: https://doi.org/10.7554/eLife.704081
 

%% Parse Calyx of Held Data

    Calyx_Calcium = DoseResp_of_Calyx(:, 1);        % Calcium Concentration for evoked peak release rate
    Calyx_PeakRelRate = DoseResp_of_Calyx(:, 2);    % Corresponding Peak Rate 
    Calyx_TTP_Calcium = TimeToPeak_of_Calyx(:, 1);  % Calcium Concentration for evoked delay
    Calyx_TTP_Time = TimeToPeak_of_Calyx(:, 2);     % Corresponding Delay Time (Measured as time to peak release rate)
  
%% Parse cMFB data 
  
    cMFB_Calcium = DoseResp_of_cMFBs(:, 1);         % Calcium Concentration for evoked peak release rate    
    cMFB_PeakRelRate = DoseResp_of_cMFBs(:, 2);     % Corresponding Peak Rate 
    cMFB_Calcium_ttp = TimeToPeak_of_cMFBs(:, 1);    % Calcium Concentration for evoked delay
    cMFB_TTP_Time = TimeToPeak_of_cMFBs(:, 2);     % Corresponding Delay Time (Measured as time to peak release rate) 
    
    
%% Simulation Variables

    Cac = sort(DoseResp_of_Calyx(:, 1));    % uM | Temporary Array of Cytosolic calcium values used as calcium step in simulation
    Can = sort(DoseResp_of_Calyx(:, 1));    % uM | Temporary Array of Nanodomain calcium values used as calcium step in simulation
    
    duration = 100;                         % ms | Duration of simulation
    Ca_stimulus_on_time = 10;               % ms | Time when step-wise calcium stimulation is applied
    tspan = linspace(0, duration, duration*100 + 1);                            % ms | simulation time-span
    tspan_1 = linspace(0, Ca_stimulus_on_time, Ca_stimulus_on_time*100 + 1);    % ms | simulation time-span before Calcium step
    tspan_2 = linspace(Ca_stimulus_on_time, duration, duration*100 + 1);        % ms | simulation time-span after calcium step 


%% Dual Sensor Model
  
% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcium Dual Sensor model parameters  |   Description

    alpha = 0.061200;     % |   /uMms   |   Association rate for synchronous release
    beta = 2.320000;      % |   /ms     |   Dissociation= rate for synchronous release
    Chi = 0.002933;       % |   /uMms   |   Association rate for Asynchronous release
    delta = 0.014829;     % |   /ms     |   Dissociation rate for Asynchronous release
    a = 0.025007;
    b_dual = 0.250007;    % |           |   Cooperativity factor
    gam_1 = 0.000009;     % |   /ms     |   Spontaneous release rate
    gam_2 = 2.000008;     % |   /ms     |   Synchronous release rate
    gam_3 = a*gam_2;      % |   /ms     |   Asynchronous release rate

    
% Reccruitment model parameters |   Unit    |   Description
    Krf_dual =  1/10.340000;  % |   /ms     |   rate of recovery of refractoriness
    Kmob_dual = 0.000050;     % |   /uMms   |   Mobilization rate
    Kdemob_dual = 0.0022;     % |   /ms     |   Demobilization rate  
    Kprime_dual = 0.027990;   % |   /uMms   |   Priming rate 
    Kupr_dual = 0.005356;     % |   /ms     |   Unpriming rate
    Kattach_dual = 0.00150;   % |   /uMms   |   Attachement rate
    Kdetach_dual = 0.001158;  % |   /ms     |   Detachhement rate
    
    
 
    DualSensorModel_Slow = cell(length(Cac), 1);            % Slow Component of release
    DualSensorModel_Fast = cell(length(Cac), 1);            % Fast Component of release
    DualSensor_SpontRate = cell(length(Cac), 1);            % Spontaneous release component
    DualSensor_PeakRecruitRate = cell(length(Cac), 1);      % Peak rate of recruitment
    DualSensorModelRate = cell(length(Cac), 1);             % Overall release rate
    DualSensorModelPeakRate = cell(length(Cac), 1);         % Peak release rate
    DualSensor_ttp = cell(length(Cac), 1);                  % Time to peak release 


%% Modified Allosteric Model

  
   % Paramter                     |     Unit    | 
    Kon = 0.097909;             % |     /uMms   | Forward reaction rate
    Koff = 3.316730;            % |     /ms     | Backward reaction rate
    I = 0.0000001;              % |     /ms     | vesicle fusion rate constant
    F = 28.693830;              % |             | Vescile fusion cooperativity open Ca2+ binding
    b_allo = 0.5008504;         % |             | Cooperativity factor
    Krf_allo = 1/6.339942;      % |     /ms     | rate of recovery of refractoriness
    Kmob_allo = 0.003858;       % |     /uMms   | Mobilization rate
    Kdemob_allo = 0.002192;     % |     /ms     | Demobilization rate 
    Kprime_allo = 0.028560;     % |     /uMms   | Priming rate 
    Kupr_allo = 0.003124;       % |     /ms     | Unpriming rate
    kattach_allo = 0.000144;    % |     /uMms   | Attachement rate
    kdetach_allo = 0.002413995; % |     /ms     | Detachhement rate
    
   
    
    AllostericModel_Slow = cell(length(Cac), 1);            % Slow Component of release           
    AllostericModel_Fast = cell(length(Cac), 1);            % Fast Component of release
    AllostericModelRate = cell(length(Cac), 1);             % Overall release rate
    AllostericModelPeakRate = cell(length(Cac), 1);         % Peak release rate
    Allosteric_PeakRecruitRate = cell(length(Cac), 1);      % Peak rate of recruitment
    Allosteric_ttp = cell(length(Cac), 1);                  % Time to peak release 



%% DoseResp Estimation

    for i=1:length(Cac) % loop though each calcium step
        
%% Allosteric Model

        % State initial values; Reserve Pool (R) - 170; Docked Pool (U) - 20; SRP (V0) - 5; FRP (W0) - 5; Others - 0                                        
        y_init_allo = [5 0 0 0 0 0 5 0 0 0 0 0 170 20 0 0];  
        
        % Solve ODEs pre-stimulation 
        [t_allo1, y_allo1] = ode15s(@(t_allo1, y_allo1)AllostericModel(y_allo1, 0, 0), tspan_1, y_init_allo);
        
        % Solve ODEs post-stimulation (Calcium step)
        [t_allo2, y_allo2] = ode15s(@(t_allo2, y_allo2)AllostericModel(y_allo2, Cac(i), Can(i)), tspan_2, y_allo1(end, :)); 
        
        y_allo = [y_allo1; y_allo2];    % concatenate pre- and post-stimulation ode solutions
        t_allo = [t_allo1; t_allo2];    % concatenate pre- and post-stimulation time steps
        
        R = y_allo(:, 13); U = y_allo(:, 14); RFv = y_allo(:, 15); RFw = y_allo(:, 16); 
        V0 = y_allo(:, 1); V1 = y_allo(:, 2); V2 = y_allo(:, 3); V3 = y_allo(:, 4); V4 = y_allo(:, 5); V5 = y_allo(:, 6);
        W0 = y_allo(:, 7); W1 = y_allo(:, 8); W2 = y_allo(:, 9); W3 = y_allo(:, 10); W4 = y_allo(:, 11); W5 = y_allo(:, 12);
        
        AllostericModel_Slow{i} = I.*(V0 + F.*V1 + (F^2).*V2 + (F^3).*V3 + (F^4).*V4 + (F^5).*V5);          

        AllostericModel_Fast{i} = I.*(W0 + F.*W1 + (F^2).*W2 + (F^3).*W3 + (F^4).*W4 + (F^5).*W5);   
                                
        AllostericModelRate{i} = AllostericModel_Slow{i} + AllostericModel_Fast{i};
        [AllostericModelPeakRate{i}, max_index_allo] = max(AllostericModelRate{i}, [], 1);
        Allosteric_ttp{i} = t_allo(max_index_allo);
        Allosteric_PeakRecruitRate{i} = max(Kmob_dual*Cac(i).*R + Kprime_dual.*U.*Cac(i).*(1 - RFv));
        

        
%% Dual Sensor Model

        % State initial values; Reserve Pool (R) - 170; Docked Pool (U) - 20; SRP (V00) - 5; FRP (W00) - 5; Others - 0                                        
            
        y_init_dual = [5 0 0 0 0 0 0 0 0 0 0 0 ...
                        0 0 0 0 0 0 5 0 0 0 0 0 0 ...
                        0 0 0 0 0 0 0 0 0 0 0 170 20 ...
                        0 0];
                    
        % Solve ODEs pre-stimulation 
        [t_dual1, y_dual1] = ode15s(@(t_dual1, y_dual1)DualSensorModel(y_dual1, 0, 0), tspan_1, y_init_dual);
        
        % Solve ODEs post-stimulation (Calcium step)        
        [t_dual2, y_dual2] = ode15s(@(t_dual2, y_dual2)DualSensorModel(y_dual2, Cac(i), Can(i)), tspan_2, y_dual1(end, :));
        
        y_dual = [y_dual1; y_dual2];    % concatenate pre- and post-stimulation ode solutions
        t_dual = [t_dual1; t_dual2];    % concatenate pre- and post-stimulation time steps
        
        
        R = y_dual(:,37) ; U = y_dual(:,38); RFv = y_dual(:,39); RFw = y_dual(:,40);

        V00 = y_dual(:,1); V01 = y_dual(:,2); V02 = y_dual(:,3); V10 = y_dual(:,4); V11 = y_dual(:,5);...
        V12 = y_dual(:,6); V20 = y_dual(:,7); V21 = y_dual(:,8); V22 = y_dual(:,9); 
        V30 = y_dual(:,10); V31 = y_dual(:,11); V32 = y_dual(:,12); V40 = y_dual(:,13); V41 = y_dual(:,14);...
        V42 = y_dual(:,15); V50 = y_dual(:,16); V51 = y_dual(:,17); V52 = y_dual(:,18);

        W00 = y_dual(:,19); W01 = y_dual(:,20); W02 = y_dual(:,21); W10 = y_dual(:,22); W11 = y_dual(:,23);...
        W12 = y_dual(:,24); W20 = y_dual(:,25); W21 = y_dual(:,26); W22 = y_dual(:,27); 
        W30 = y_dual(:,28); W31 = y_dual(:,29); W32 = y_dual(:,30); W40 = y_dual(:,31); W41 = y_dual(:,32);...
        W42 = y_dual(:,33); W50 = y_dual(:,34); W51 = y_dual(:,35); W52 = y_dual(:,36);

        DualSensorModel_Slow{i}= gam_1.*V00 + gam_2.*(V50 + V51 + V52) + gam_3.*(V02 + V12 + V22 + V32 + V42 + V52);          

        DualSensorModel_Fast{i} = gam_1.*W00 + gam_2.*(W50 + W51 + W52) + gam_3.*(W02 + W12 + W22 + W32 + W42 + W52);  
        
        DualSensorModelRate{i} = DualSensorModel_Slow{i} + DualSensorModel_Fast{i};
        [DualSensorModelPeakRate{i}, max_index_dual] = max(DualSensorModelRate{i}, [], 1);
        DualSensor_ttp{i} = t_dual(max_index_dual);
        DualSensor_PeakRecruitRate{i} = max(Kmob_allo*Cac(i).*R + Kprime_allo.*U.*Cac(i).*(1 - RFv));
        DualSensor_SpontRate{i} = gam_1.*V00 + gam_1.*W00;
        
    end 


    
%% Fit peak release rate from both models using Hill Equations located at the bottom of file


%% Alloesteric Dose Response(HILL) Fit

params_init_allo = [-0.002934 6.724 0.9638 4.027];  % Initial paramters
lb_allo = [];                                       % Lower bound for paramters if required
ub_allo = [];                                       % Upper bound for paramters if required

x = log10(Cac);
Y_allo = cell2mat(AllostericModelPeakRate);
[params_allo] = lsqcurvefit(@(params_allo, x)DoseResponseFit(params_allo, x),params_init_allo, x, Y_allo, lb_allo,ub_allo);
                                      
AllostericPeakRate_fit = DoseResponseFit(params_allo, log10(Cac));


%% Dualsensor Dose Response(HILL) Fit


params_init_dual = [-0.01724 6.401 0.9866 2.524];   % Initial paramters
lb_dual = [];                                       % Lower bound for paramters if required
ub_dual = [];                                       % Upper bound for paramters if required 

Y_dual = cell2mat(DualSensorModelPeakRate);
[params_dual] = lsqcurvefit(@(params_dual, x)DoseResponseFit(params_dual, x),params_init_dual, x, Y_dual, lb_dual,ub_dual);
                                      
DualSensorPeakRate_fit = DoseResponseFit(params_dual, log10(Cac));


%% Visualize Spontaneous Release Rate for Resting Level Calcium

SpontRate = SpontaneousRelRate_from_Suhita_2010(:, 2);
time_SpontRate = SpontaneousRelRate_from_Suhita_2010(:, 1);
index = find(round(Cac, 2) == 0.1, 1);
figure
%semilogy(time_SpontRate, SpontRate , "k.", "Color",'#A2142F','LineWidth', 1, 'MarkerSize', 10)
%hold on
semilogy(t_dual, DualSensor_SpontRate{index}, "-", "Color",'#D95319', 'LineWidth', 1, 'MarkerSize', 10)
hold on
legend({'Dual Sensor (sim)'},'Location', 'northwest', 'FontSize',6)
ylabel('Release Rate (ms^{-1})','FontSize',6,'FontWeight','bold','Color','k')
xlabel('time (ms)','FontSize',6,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
hold off


%Get Current Figure (GCF) & Set image size before saving image
width = 7;  % cm 
height = 6; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 SpontaneousRate_DataAndSimulation_vs_Time
saveas(gcf, 'SpontaneousRate_DataAndSimulation_vs_Time', 'fig')

hold off


movefile SpontaneousRate_DataAndSimulation_vs_Time.png ../../results/Figure_1
movefile SpontaneousRate_DataAndSimulation_vs_Time.fig ../../results/Figure_1

%% Visualize Release Rate for Each Calcium Step 

%%%%%%%%%%%%%%%%%%%%%%%%   Allosteric Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

pos = [0.1 0.4 0.7 0.55];
subplot("Position", pos)
ColorSet = jet(length(Cac));

hold all;
for i=1:length(Cac)
    plot(t_allo, AllostericModelRate{i}, 'LineWidth', 1, 'MarkerSize', 10)
end

%title('Total release from Allosteric model following Ca^{2+} clamp','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Release Rate (vesicles/ms)','FontSize',6,'FontWeight','bold','Color','k')
xlabel('time (ms)','FontSize',6,'FontWeight','bold','Color','k')
set(gca, 'ColorOrder', ColorSet);
set(gca, 'Colormap', ColorSet);
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colorbar('Ticks',[0.5, 6, 12],...
         'TickLabels',{'0.5 (\muM)','6','12  (\muM)'})
caxis([min(Cac) max(Cac)])
set(gca, 'box', 'off')
hold off


pos = [0.1 0.1 0.575 0.15];
subplot("Position", pos)
ColorSet = jet(length(Cac));
Clamped_Calcium = zeros(length(t_allo), 1);
set(gca, 'ColorOrder', ColorSet);
hold all;
for i=1:length(Cac)
    
    for j=1:length(t_allo)
        if t_allo(j) < Ca_stimulus_on_time
            Clamped_Calcium(j) = 0;
        elseif t_allo(j) >= Ca_stimulus_on_time
            Clamped_Calcium(j) = Cac(i);
        end
    end
    
    plot(t_allo, Clamped_Calcium, 'LineWidth', 1, 'MarkerSize', 10)
end

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
str = {'Ca^{2+} steps'};
text(110, 6, str, 'FontSize',5,'Color','k')
set(gca, 'box', 'off')
hold off

%Get Current Figure (GCF) & Set image size before saving image
width = 7;  % cm 
height = 6; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Clamped_AlloReleaseRate_vs_Time
saveas(gcf, 'Clamped_AlloReleaseRate_vs_Time', 'fig')
hold off

movefile Clamped_AlloReleaseRate_vs_Time.png ../../results/Figure_1
movefile Clamped_AlloReleaseRate_vs_Time.fig ../../results/Figure_1

%%%%%%%%%%%%%%%%%%%%%%%%  Dual-Sensor Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

pos = [0.1 0.4 0.7 0.55];
subplot("Position", pos)
ColorSet = jet(length(Cac));

hold all;
for i=1:length(Cac)
    plot(t_allo, DualSensorModelRate{i}, 'LineWidth', 1, 'MarkerSize', 10)
end

%title('Total release from Dual-sensor model following Ca^{2+} clamp','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Release Rate (vesicles/ms)','FontSize',6,'FontWeight','bold','Color','k')
xlabel('time (ms)','FontSize',6,'FontWeight','bold','Color','k')
set(gca, 'ColorOrder', ColorSet);
set(gca, 'Colormap', ColorSet);
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colorbar('Ticks',[0.5, 6, 12],...
         'TickLabels',{'0.5 (\muM)','6','12  (\muM)'})
caxis([min(Cac) max(Cac)])
set(gca, 'box', 'off')
hold off


pos = [0.1 0.1 0.575 0.15];
subplot("Position", pos)
ColorSet = jet(length(Cac));
Clamped_Calcium = zeros(length(t_allo), 1);
set(gca, 'ColorOrder', ColorSet);
hold all;

for i=1:length(Cac)
    
    for j=1:length(t_allo)
        if t_allo(j) < Ca_stimulus_on_time
            Clamped_Calcium(j) = 0;
        elseif t_allo(j) >= Ca_stimulus_on_time
            Clamped_Calcium(j) = Cac(i);
        end
    end
    
    plot(t_allo, Clamped_Calcium, 'LineWidth', 1, 'MarkerSize', 10)
end

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
str = {'Ca^{2+} steps'};
text(110, 6, str, 'FontSize',5,'Color','k')
set(gca, 'box', 'off')
hold off

%Get Current Figure (GCF) & Set image size before saving image
width = 7;  % cm 
height = 6; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Clamped_DualReleaseRate_vs_Time
saveas(gcf, 'Clamped_DualReleaseRate_vs_Time', 'fig')
hold off

movefile Clamped_DualReleaseRate_vs_Time.png ../../results/Figure_1
movefile Clamped_DualReleaseRate_vs_Time.fig ../../results/Figure_1


%% Time To Peak Rate (ms)

Ca_temp = Cac(Cac > 2);                         % Selection criteria required for filtering Ca steps with approx no delay. 
Allosteric_ttp2 = cell2mat(Allosteric_ttp);
Allo_ttp_temp = Allosteric_ttp2(Cac > 2);

DualSensor_ttp2 = cell2mat(DualSensor_ttp);
Dual_ttp_temp = DualSensor_ttp2(Cac > 2); 
    
    
figure

plot(Calyx_TTP_Calcium, Calyx_TTP_Time ,".", "Color",'#A2142F','LineWidth', 1, 'MarkerSize', 10)
hold on
plot(Ca_temp(Allo_ttp_temp < 60), Allo_ttp_temp(Allo_ttp_temp < 60) - Ca_stimulus_on_time,...
     "-", "Color",'#D95319', 'LineWidth', 1, 'MarkerSize', 10)
hold on
plot(Ca_temp(Dual_ttp_temp < 60), Dual_ttp_temp(Dual_ttp_temp < 60) - Ca_stimulus_on_time,...
    "-", "Color",'#0072BD', 'LineWidth', 1, 'MarkerSize', 10)
hold on
plot(cMFB_Calcium_ttp, cMFB_TTP_Time ,".", "Color", '#000000', 'LineWidth', 1, 'MarkerSize', 10)
hold on

legend({"Calyx of Held", 'Allosteric', 'Dual Sensor'},'Location', 'northeast', 'FontSize',6)
ylabel('Time to peak release (ms)','FontSize',6,'FontWeight','bold','Color','k')
xlabel('Ca^{2+} (\muM)','FontSize',6,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')

hold off


%Get Current Figure (GCF) & Set image size before saving image
width = 7;  % cm 
height = 6; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng  -r1000 TimeToPeakRelease_vs_Calcium
saveas(gcf, 'TimeToPeakRelease_vs_Calcium', 'fig')
hold off


movefile TimeToPeakRelease_vs_Calcium.png ../../results/Figure_1
movefile TimeToPeakRelease_vs_Calcium.fig ../../results/Figure_1

%% Dose Response Curve


cMFB_Calcium = DoseResp_of_cMFBs(:, 1); 
figure
loglog(Calyx_Calcium , Calyx_PeakRelRate, ".", "Color",'#A2142F','LineWidth', 1, 'MarkerSize', 5)
hold on
loglog(Cac, cell2mat(AllostericModelPeakRate), "-", "Color",'#D95319', 'LineWidth', 1, 'MarkerSize', 10)
hold on
loglog(Cac, cell2mat(DualSensorModelPeakRate), "-", "Color",'#0072BD', 'LineWidth', 1, 'MarkerSize', 10)
hold on
loglog(cMFB_Calcium, cMFB_PeakRelRate, "^", "Color", '#000000', 'LineWidth', 1, 'MarkerSize', 1)
hold on
str = {'Dual Sensor:',...
        strcat('EC_{50} =', " ", num2str(round(10^(params_dual(3)))), " ", "(\muM)"),...
        strcat('Hill Slope =', " ", num2str(round(params_dual(4))))};
text(0.02, 7, str, 'FontSize',5,'Color','k')
legend({"Calyx of Held", 'Allosteric', 'Dual Sensor', 'cMFBs | Eshra et al.  2021'},'Location', 'northwest', 'FontSize',3)
ylabel('Peak release rate (vesicle ms^{-1})','FontSize',6,'FontWeight','bold','Color','k')
xlabel('Ca^{2+} (\muM)','FontSize',6,'FontWeight','bold','Color','k')
set(gca, 'box', 'off')
hold off


%Get Current Figure (GCF) & Set image size before saving image
width = 7;  % cm 
height = 6; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 PeakRelease_vs_Calcium
saveas(gcf, 'PeakRelease_vs_Calcium', 'fig')
hold off


movefile PeakRelease_vs_Calcium.png ../../results/Figure_1
movefile PeakRelease_vs_Calcium.fig ../../results/Figure_1

 %% Release ODEs
 
 
%%%%%%%%%%%%%%%%%%%%%%%%   Allosteric Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = AllostericModel(y, Cac, Ca_vgcc)

%{
Allosteric modulation of the presynaptic Ca21 sensor for vesicle fusion
Xuelin Lou1, Volker Scheuss1? & Ralf Schneggenburger1
Calyx Of Held
%}

   % Paramter                     |     Unit    | 
    Kon = 0.097909;             % |     /uMms   | Forward reaction rate
    Koff = 3.316730;            % |     /ms     | Backward reaction rate
    I = 0.0000001;              % |     /ms     | vesicle fusion rate constant
    F = 28.693830;              % |             | Vescile fusion cooperativity open Ca2+ binding
    b_allo = 0.5008504;         % |             | Cooperativity factor
    Krf_allo = 1/6.339942;      % |     /ms     | rate of recovery of refractoriness
    Kmob_allo = 0.003858;       % |     /uMms   | Mobilization rate
    Kdemob_allo = 0.002192;     % |     /ms     | Demobilization rate 
    Kprime_allo = 0.028560;     % |     /uMms   | Priming rate 
    Kupr_allo = 0.003124;       % |     /ms     | Unpriming rate
    kattach_allo = 0.000144;    % |     /uMms   | Attachement rate
    kdetach_allo = 0.002413995; % |     /ms     | Detachhement rate
    

    V0 = y(1);
    V1 = y(2);
    V2 = y(3);
    V3 = y(4);
    V4 = y(5);
    V5 = y(6); 
     
    W0 = y(7);
    W1 = y(8);
    W2 = y(9);
    W3 = y(10);
    W4 = y(11);
    W5 = y(12);
     
    R_allo = y(13);
    U_allo = y(14);
    RFv_allo = y(15);
    RFw_allo = y(16); 

    Vtotal_allo = V0 + V1 + V2 + V3 + V4 + V5; % Total vesicles in SRP

    Wtotal_allo = W0 + W1 + W2 + W3 + W4 + W5; % Total vesicles in FRP

    dR_allodt = -Kmob_allo*Cac*R_allo + U_allo*Kdemob_allo;
    dU_allodt =  -Kdemob_allo*U_allo + Kmob_allo*Cac*R_allo - Kprime_allo*U_allo*Cac*(1 - RFv_allo) + Kupr_allo*V0;

    dRFv_allodt = (1 - RFv_allo)*(I*(V0 + V1*F + V2*F^2 + V3*F^3 + V4*F^4 + V5*F^5))/Vtotal_allo - Krf_allo*RFv_allo;

    dRFw_allodt = (1 - RFw_allo)*(I*(W0 + W1*F + W2*F^2 + W3*F^3 + W4*F^4 + W5*F^5))/Wtotal_allo - Krf_allo*RFw_allo;

    dV0dt = Kprime_allo*U_allo*Cac*(1 - RFv_allo) - Kupr_allo*V0 - kattach_allo*V0*Ca_vgcc*(1 - RFw_allo) +...
            kdetach_allo*W0 + Koff*V1 - I*V0 - 5*Kon*Cac*V0; 
    dV1dt = 2*Koff*b_allo*V2 - Koff*V1 + 5*Kon*Cac*V0 - 4*Kon*Cac*V1 - I*F*V1;
    dV2dt = 3*Koff*V3*b_allo^2 - 2*Koff*V2*b_allo + 4*Kon*Cac*V1 - 3*Kon*Cac*V2 - I*V2*F^2;
    dV3dt = 4*Koff*V4*b_allo^3 - 3*Koff*V3*b_allo^2 + 3*Kon*Cac*V2 - 2*Kon*Cac*V3 - I*V3*F^3;
    dV4dt = 5*Koff*V5*b_allo^4 - 4*Koff*V4*b_allo^3 + 2*Kon*Cac*V3 - Kon*Cac*V4 - I*V4*F^4;
    dV5dt = Kon*Cac*V4 - 5*Koff*V5*b_allo^4 - I*V5*F^5;


    dW0dt = kattach_allo*W0*Ca_vgcc*(1 - RFw_allo) - kdetach_allo*W0 + Koff*W1 - I*W0 - 5*Kon*Ca_vgcc*W0; 
    dW1dt = 2*Koff*b_allo*W2 - Koff*W1 + 5*Kon*Ca_vgcc*W0 - 4*Kon*Ca_vgcc*W1 - I*F*W1;
    dW2dt = 3*Koff*W3*b_allo^2 - 2*Koff*W2*b_allo + 4*Kon*Ca_vgcc*W1 - 3*Kon*Ca_vgcc*W2 - I*W2*F^2;
    dW3dt = 4*Koff*W4*b_allo^3 - 3*Koff*W3*b_allo^2 + 3*Kon*Ca_vgcc*W2 - 2*Kon*Ca_vgcc*W3 - I*W3*F^3;
    dW4dt = 5*Koff*W5*b_allo^4 - 4*Koff*W4*b_allo^3 + 2*Kon*Ca_vgcc*W3 - Kon*Ca_vgcc*W4 - I*W4*F^4;
    dW5dt = Kon*Ca_vgcc*W4 - 5*Koff*W5*b_allo^4 - I*W5*F^5;
    
    
    dydt=[dV0dt; dV1dt; dV2dt; dV3dt; dV4dt; dV5dt; dW0dt; dW1dt; dW2dt; dW3dt; dW4dt; dW5dt;...
          dR_allodt; dU_allodt; dRFv_allodt; dRFw_allodt];

end


%%%%%%%%%%%%%%%%%%%%%%%%  Dual-Sensor Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dydt_5 = DualSensorModel(y, Cac, Ca_vgcc)

%{
A dual sensor model for neurotransmitter release in central synapse
Sun et al. 2007
Calyx Of Held
%}

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcium Dual Sensor model parameters  |   Description

    alpha = 0.061200;     % |   /uMms   |   Association rate for synchronous release
    beta = 2.320000;      % |   /ms     |   Dissociation= rate for synchronous release
    Chi = 0.002933;       % |   /uMms   |   Association rate for Asynchronous release
    delta = 0.014829;     % |   /ms     |   Dissociation rate for Asynchronous release
    a = 0.025007;
    b_dual = 0.250007;    % |           |   Cooperativity factor
    gam_1 = 0.000009;     % |   /ms     |   Spontaneous release rate
    gam_2 = 2.000008;     % |   /ms     |   Synchronous release rate
    gam_3 = a*gam_2;      % |   /ms     |   Asynchronous release rate

    
% Reccruitment model parameters |   Unit    |   Description
    Krf_dual =  1/10.340000;  % |   /ms     |   rate of recovery of refractoriness
    Kmob_dual = 0.000050;     % |   /uMms   |   Mobilization rate
    Kdemob_dual = 0.0022;     % |   /ms     |   Demobilization rate  
    Kprime_dual = 0.027990;   % |   /uMms   |   Priming rate 
    Kupr_dual = 0.005356;     % |   /ms     |   Unpriming rate
    Kattach_dual = 0.00150;   % |   /uMms   |   Attachement rate
    Kdetach_dual = 0.001158;  % |   /ms     |   Detachhement rate
    
    
    V00 = y(1); 
    V01 = y(2);
    V02 = y(3);
    V10 = y(4);
    V11 = y(5);
    V12 = y(6);
    V20 = y(7);
    V21 = y(8);
    V22 = y(9);
    V30 = y(10);
    V31 = y(11);
    V32 = y(12);
    V40 = y(13);
    V41 = y(14);
    V42 = y(15);
    V50 = y(16);
    V51 = y(17);
    V52 = y(18);
    
    W00 = y(19); 
    W01 = y(20);
    W02 = y(21);
    W10 = y(22);
    W11 = y(23);
    W12 = y(24);
    W20 = y(25);
    W21 = y(26);
    W22 = y(27);
    W30 = y(28);
    W31 = y(29);
    W32 = y(30);
    W40 = y(31);
    W41 = y(32);
    W42 = y(33);
    W50 = y(34);
    W51 = y(35);
    W52 = y(36);
    
    R_dual = y(37);
    U_dual = y(38);
    RFv_dual = y(39);
    RFw_dual= y(40); 
    
    Vtotal_dual = V00 + V01 + V02 + ...
              V10 + V11 + V12 + ...
              V20 + V21 + V22 + ...         % Total vesicles in the SRP
              V30 + V31 + V32 + ...
              V40 + V41 + V42 + ...
              V50 + V51 + V52; 

    Wtotal_dual = W00 + W01 + W02 + ...
              W10 + W11 + W12 + ...
              W20 + W21 + W22 + ...         % Total vesicles in the FRP
              W30 + W31 + W32 + ...
              W40 + W41 + W42 + ...
              W50 + W51 + W52;   

    dR_dualdt = -Kmob_dual*Cac*R_dual + U_dual*Kdemob_dual;
    dU_dualdt = -Kprime_dual*U_dual*Cac*(1 - RFv_dual) + Kupr_dual*V00;

    dRFv_dualdt = ((1 - RFv_dual)*(gam_3*(V02 + V12 + V22 + V32 + V42 + V52) + ...
                        gam_2*(V50 + V51 + V52) + gam_1*V00)/Vtotal_dual - Krf_dual*RFv_dual);

    dRFw_dualdt = ((1 - RFw_dual)*(gam_3*(W02 + W12 + W22 + W32 + W42 + W52) + ...
                        gam_2*(W50 + W51 + W52) + gam_1*W00)/Wtotal_dual - Krf_dual*RFw_dual);

 
    dV00 = Kprime_dual*U_dual*Cac*(1 - RFv_dual) - Kupr_dual*V00 - Kattach_dual*V00*Ca_vgcc*(1 - RFw_dual) + Kdetach_dual*W00 + ...
             beta*V10 - 5*alpha*V00*Cac + delta*V01 - 2*Chi*V00*Cac - gam_1*V00;
    dV01 = 2*Chi*V00*Cac - delta*V01 + 2*b_dual*delta*V02 - Chi*Cac*V01 + beta*V11 - 5*alpha*Cac*V01;
    dV02 = Chi*Cac*V01 - 2*b_dual*delta*V02 - gam_3*V02 - 5*alpha*Cac*V02 + beta*V12;
    dV10 = 5*alpha*Cac*V00 - beta*V10 - 4*alpha*Cac*V10 + 2*b_dual*beta*V20 - 2*Chi*Cac*V10 + delta*V11;
    dV11 = 5*alpha*Cac*V01 - beta*V11 - 4*alpha*Cac*V11 + 2*b_dual*beta*V21 + 2*Chi*Cac*V10 - delta*V11 - Chi*Cac*V11 + 2*b_dual*delta*V12;
    dV12 = 5*alpha*Cac*V02 - beta*V12 - 4*alpha*Cac*V12 + 2*b_dual*beta*V22 + Chi*Cac*V11 - 2*b_dual*delta*V12 - gam_3*V12;
    dV20 = 4*alpha*Cac*V10 - 2*b_dual*beta*V20 - 3*alpha*Cac*V20 + 3*(b_dual^2)*beta*V30 - 2*Chi*Cac*V20 + delta*V21;
    dV21 = 4*alpha*Cac*V11 - 2*b_dual*beta*V21 - 3*alpha*Cac*V21 + 3*(b_dual^2)*beta*V31 + 2*Chi*Cac*V20 - delta*V21 - Chi*Cac*V21 + 2*b_dual*delta*V22;
    dV22 = 4*alpha*Cac*V12 - 2*b_dual*beta*V22 - 3*alpha*Cac*V22 + 3*(b_dual^2)*beta*V32 + Chi*Cac*V21 - 2*b_dual*delta*V22 - gam_3*V22;
    dV30 = 3*alpha*Cac*V20 - 3*(b_dual^2)*beta*V30 - 2*alpha*Cac*V30 + 4*(b_dual^3)*beta*V40 - 2*Chi*Cac*V30 + delta*V31;
    dV31 = 3*alpha*Cac*V21 - 3*(b_dual^2)*beta*V31 - 2*alpha*Cac*V31 + 4*(b_dual^3)*beta*V41 + 2*Chi*Cac*V30 - delta*V31 - Chi*Cac*V31 + 2*b_dual*delta*V32;
    dV32 = 3*alpha*Cac*V22 - 3*(b_dual^2)*beta*V32 - 2*alpha*Cac*V32 + 4*(b_dual^2)*beta*V42 + Chi*Cac*V31 - 2*b_dual*delta*V32 - gam_3*V32;
    dV40 = 2*alpha*Cac*V30 - 4*(b_dual^3)*beta*V40 - alpha*Cac*V40 + 5*(b_dual^4)*beta*V50 - 2*Chi*Cac*V40 + delta*V41;
    dV41 = 2*alpha*Cac*V31 - 4*(b_dual^3)*beta*V41 - alpha*Cac*V41 + 5*(b_dual^4)*beta*V51 + 2*Chi*Cac*V40 - delta*V41 - Chi*Cac*V41 + 2*b_dual*delta*V42;
    dV42 = 2*alpha*Cac*V32 - 4*(b_dual^3)*beta*V42 - alpha*Cac*V42 + 5*(b_dual^4)*beta*V52 + Chi*Cac*V41 - 2*b_dual*delta*V42 - gam_3*V42;
    dV50 = alpha*Cac*V40 - 5*(b_dual^4)*beta*V50 - 2*Chi*Cac*V50 + delta*V51 - gam_2*V50;
    dV51 = alpha*Cac*V41 - 5*(b_dual^4)*beta*V51 + 2*Chi*Cac*V50 - delta*V51 - Chi*Cac*V51 + 2*b_dual*delta*V52 - gam_2*V51;
    dV52 = alpha*Cac*V42 - 5*(b_dual^4)*beta*V52 + Chi*Cac*V51 - 2*b_dual*delta*V52 - gam_3*V52 - gam_2*V52;


    dW00 = Kattach_dual*V00*Ca_vgcc*(1 - RFw_dual) - Kdetach_dual*W00 + ...
             beta*W10 - 5*alpha*W00*Ca_vgcc + delta*W01 - 2*Chi*W00*Ca_vgcc - gam_1*W00;
    dW01 = 2*Chi*W00*Ca_vgcc - delta*W01 + 2*b_dual*delta*W02 - Chi*Ca_vgcc*W01 + beta*W11 - 5*alpha*Ca_vgcc*W01;
    dW02 = Chi*Ca_vgcc*W01 - 2*b_dual*delta*W02 - gam_3*W02 - 5*alpha*Ca_vgcc*W02 + beta*W12;
    dW10 = 5*alpha*Ca_vgcc*W00 - beta*W10 - 4*alpha*Ca_vgcc*W10 + 2*b_dual*beta*W20 - 2*Chi*Ca_vgcc*W10 + delta*W11;
    dW11 = 5*alpha*Ca_vgcc*W01 - beta*W11 - 4*alpha*Ca_vgcc*W11 + 2*b_dual*beta*W21 + 2*Chi*Ca_vgcc*W10 - delta*W11 - Chi*Ca_vgcc*W11 + 2*b_dual*delta*W12;
    dW12 = 5*alpha*Ca_vgcc*W02 - beta*W12 - 4*alpha*Ca_vgcc*W12 + 2*b_dual*beta*W22 + Chi*Ca_vgcc*W11 - 2*b_dual*delta*W12 - gam_3*W12;
    dW20 = 4*alpha*Ca_vgcc*W10 - 2*b_dual*beta*W20 - 3*alpha*Ca_vgcc*W20 + 3*(b_dual^2)*beta*W30 - 2*Chi*Ca_vgcc*W20 + delta*W21;
    dW21 = 4*alpha*Ca_vgcc*W11 - 2*b_dual*beta*W21 - 3*alpha*Ca_vgcc*W21 + 3*(b_dual^2)*beta*W31 + 2*Chi*Ca_vgcc*W20 - delta*W21 - Chi*Ca_vgcc*W21 + 2*b_dual*delta*W22;
    dW22 = 4*alpha*Ca_vgcc*W12 - 2*b_dual*beta*W22 - 3*alpha*Ca_vgcc*W22 + 3*(b_dual^2)*beta*W32 + Chi*Ca_vgcc*W21 - 2*b_dual*delta*W22 - gam_3*W22;
    dW30 = 3*alpha*Ca_vgcc*W20 - 3*(b_dual^2)*beta*W30 - 2*alpha*Ca_vgcc*W30 + 4*(b_dual^3)*beta*W40 - 2*Chi*Ca_vgcc*W30 + delta*W31;
    dW31 = 3*alpha*Ca_vgcc*W21 - 3*(b_dual^2)*beta*W31 - 2*alpha*Ca_vgcc*W31 + 4*(b_dual^3)*beta*W41 + 2*Chi*Ca_vgcc*W30 - delta*W31 - Chi*Ca_vgcc*W31 + 2*b_dual*delta*W32;
    dW32 = 3*alpha*Ca_vgcc*W22 - 3*(b_dual^2)*beta*W32 - 2*alpha*Ca_vgcc*W32 + 4*(b_dual^2)*beta*W42 + Chi*Ca_vgcc*W31 - 2*b_dual*delta*W32 - gam_3*W32;
    dW40 = 2*alpha*Ca_vgcc*W30 - 4*(b_dual^3)*beta*W40 - alpha*Ca_vgcc*W40 + 5*(b_dual^4)*beta*W50 - 2*Chi*Ca_vgcc*W40 + delta*W41;
    dW41 = 2*alpha*Ca_vgcc*W31 - 4*(b_dual^3)*beta*W41 - alpha*Ca_vgcc*W41 + 5*(b_dual^4)*beta*W51 + 2*Chi*Ca_vgcc*W40 - delta*W41 - Chi*Ca_vgcc*W41 + 2*b_dual*delta*W42;
    dW42 = 2*alpha*Ca_vgcc*W32 - 4*(b_dual^3)*beta*W42 - alpha*Ca_vgcc*W42 + 5*(b_dual^4)*beta*W52 + Chi*Ca_vgcc*W41 - 2*b_dual*delta*W42 - gam_3*W42;
    dW50 = alpha*Ca_vgcc*W40 - 5*(b_dual^4)*beta*W50 - 2*Chi*Ca_vgcc*W50 + delta*W51 - gam_2*W50;
    dW51 = alpha*Ca_vgcc*W41 - 5*(b_dual^4)*beta*W51 + 2*Chi*Ca_vgcc*W50 - delta*W51 - Chi*Ca_vgcc*W51 + 2*b_dual*delta*W52 - gam_2*W51;
    dW52 = alpha*Ca_vgcc*W42 - 5*(b_dual^4)*beta*W52 + Chi*Ca_vgcc*W51 - 2*b_dual*delta*W52 - gam_3*W52 - gam_2*W52;

    
    dydt_5 = [dV00; dV01; dV02; dV10; dV11; dV12; dV20; dV21; dV22; dV30; dV31; dV32;...
                dV40; dV41; dV42; dV50; dV51; dV52; dW00; dW01; dW02; dW10; dW11; dW12; dW20;...
                dW21; dW22; dW30; dW31; dW32; dW40; dW41; dW42; dW50; dW51; dW52;...
                dR_dualdt; dU_dualdt; dRFv_dualdt; dRFw_dualdt];
            
            
end

%%%%%%%%%%%%%%%%%%%%%%%%  Hill equation for fitting release rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y_real = DoseResponseFit(a, x)
        
        %% Gaussian Fit
        yD = a(1) + (a(2) - a(1))./(1 + 10.^((a(3) - x).*a(4)));
        y_real= yD;   
end