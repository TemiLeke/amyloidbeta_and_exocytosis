function [Po_ip3r, IcaPQ] = kinetic_schemes(y, cell_condition)
global dt N_ip3r N_vgcc state_ip3r state_PQ 

Ca_ip3 = y(2);          %               Nanodomain Calcium
IP3 = y(9);             % uM            IP3 Concentration   
Cout = 2000;            % uM            Extracellular Ca concentration, unit: uM
Ca_vgcc = y(11);        % uM            VGCC Nanodomain concentration 
V = y(4);               % mV            Membrane potential
F = 96485.33;           % Col/mole      Faraday's constant, unit: col/mole
T = 300;                % K             Temperature, unit:Kelvin
Ro = 8.31;              % J/(mole*K)    Ideal gas constant, unit: J/(mole*K)
z = 2;                  %               valence of Ca ion



%% Nernst potential for Ca %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eca = (Ro*T/(z*F))*1000*log(Cout/Ca_vgcc);           %mV 55;  

%% IP3 KINETIC SCHEME

%% WT model paramters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cell_condition=="WT"
    
    a1 = 17.050543;             % /(uM)^2
    nO = 2.473407;                           
    Kod = 0.909078;             % (uM)
    a2 = 18.49186;              % /(uM)^2
    nA = 0.093452;                        
    Kad = 1.955650e03;          % (uM)
    a3 = 2.73028e02;            % /(uM)^5
    nI = 56.84823;                      
    Kid = 0.089938;             % (uM)
    j01 = 3.031635e02;          % /(uM)^5*ms
    j12 = 3.230063e02;          % /((uM)^2)*ms  
    j22 = 4.814111;             % /((uM)^2)*ms  
    j23 = 5.356155;             % /((uM)^3)*ms  
    j45 = 5.625616;             % /((uM)^5)*ms 
    J01_tilda = 3.013284e02;    % /(uM)*ms 
    J45_tilda = 2.648741;       % /((uM)^5)*ms 
    
    
%% AD Model Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif cell_condition == "AD"
   
    a1 = 1.108278e02;           % /(uM)^2
    nO = 2.473407;                           
    Kod = 0.909078;             % (uM)
    a2 = 18.49186;              % /(uM)^2
    nA = 0.093452;                        
    Kad = 1.955650e03;          % (uM)
    a3 = 1.4041556e02;          % /(uM)^5
    nI = 56.84823;                      
    Kid = 0.089938;             % (uM)
    j01 = 3.031635e02;          % /(uM)^5*ms
    j12 = 3.230063e02;          % /((uM)^2)*ms  
    j22 = 5.3978052;            % /((uM)^2)*ms  
    j23 = 2.0652269e03;         % /((uM)^3)*ms  
    j45 = 5.4319289;            % /((uM)^5)*ms 
    J01_tilda = 3.013284e02;    % /(uM)*ms 
    J45_tilda = 8.512829e-08;   % /((uM)^5)*ms 
    
end


%% Transition Rates and Open Probability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
KO = (a1*IP3^nO)/(IP3^nO + Kod^nO);
KI = (a3*IP3^nI)/(IP3^nI + Kid^nI);
KA = (a2*IP3^nA)/(IP3^nA + Kad^nA);

kRA = ((1/(j01*Ca_ip3)) + (1/(j12*(Ca_ip3^2))))^-1;
kAR = (KA*(Ca_ip3^2)*((1/(j01*Ca_ip3)) + (1/(j12*(Ca_ip3^2)))))^-1; 

kAO = (KA*(Ca_ip3^2)*(1/(j22*(Ca_ip3^2))))^-1;
kOA = (KO*(Ca_ip3^2)*(1/(j22*(Ca_ip3^2))))^-1;

kOI = ((KO*Ca_ip3^2)*(1/(j23*Ca_ip3^3) + 1/(j45*Ca_ip3^5)))^-1;
kIO = ((KI*Ca_ip3^5)*(1/(j23*Ca_ip3^3) + 1/(j45*Ca_ip3^5)))^-1;

kRI = (1/(J01_tilda*Ca_ip3) + 1/(J45_tilda*Ca_ip3^5))^-1;
kIR = ((KI*Ca_ip3^5)*(1/(J01_tilda*Ca_ip3) + 1/(J45_tilda*Ca_ip3^5)))^-1;


%% For Stochastic IP3R Channel Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pRA = kRA*dt;
    pAR = kAR*dt;
    pOA = kOA*dt;
    pAO = kAO*dt;
    pRI = kRI*dt;
    pIR = kIR*dt;
    pIO = kIO*dt;
    pOI = kOI*dt;
    
% R = 1, A = 2, O = 3, I = 4

    OpenNPR = 0;
    for i = 1:N_ip3r
        z = rand;
        if (state_ip3r(i) == 1)
            if(z < pRA)
                state_ip3r(i) = 2;
            elseif(z >= pRA && z <  (pRA+pRI))
                state_ip3r(i) = 4;
            end
        elseif (state_ip3r(i) == 2)
            if(z < pAR)
                state_ip3r(i) = 1;
            elseif(z >= pAR && z <  (pAR+pAO))
                state_ip3r(i) = 3;
            end
        elseif (state_ip3r(i) == 3)
            if(z < pOA)
                state_ip3r(i) = 2;
            elseif( z >= pOA && z <  (pOA+pOI))
                state_ip3r(i) = 4;
            end
        elseif (state_ip3r(i) == 4)
            if(z < pIR)
                state_ip3r(i) = 1;
            elseif( z >= pIR && z <  (pIR+pIO))
                state_ip3r(i) = 3;
            end
        end
        if(state_ip3r(i) == 3)
            OpenNPR = OpenNPR + 1;
        end
    end

    Po_ip3r = (OpenNPR / N_ip3r);




%%  Stochastic VGCC channel simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Suhita Nadkarni P-Q- Type current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gPQ = 2.7e-09;                 % mS
    
    alpha_10PQ = 4.04;             % /ms 
    beta_10PQ = 2.88;              % /ms
    k1_PQ = 49.14;                 % mV 
    alpha_20PQ = 6.70;             % /ms 
    beta_20PQ = 6.30;              % /ms
    k2_PQ = 42.08;                 % mV 
    alpha_30PQ = 4.39;             % /ms 
    beta_30PQ = 8.16;              % /ms     
    k3_PQ = 55.31;                 % mV  
    alpha_40PQ = 17.33;            % /ms 
    beta_40PQ = 1.84;              % /ms    
    k4_PQ = 26.55;                 % mV
   
    alpha_1PQ = alpha_10PQ*exp(V/k1_PQ);                         % /ms 
    beta_1PQ = beta_10PQ*exp(-V/k1_PQ);                          % /ms 
    alpha_2PQ = alpha_20PQ*exp(V/k2_PQ);                         % /ms 
    beta_2PQ = beta_20PQ*exp(-V/k2_PQ);                          % /ms 
    alpha_3PQ = alpha_30PQ*exp(V/k3_PQ);                         % /ms 
    beta_3PQ = beta_30PQ*exp(-V/k3_PQ);                          % /ms 
    alpha_4PQ = alpha_40PQ*exp(V/k4_PQ);                         % /ms 
    beta_4PQ = beta_40PQ*exp(-V/k4_PQ);                          % /ms 
    
    OpenPQ = 0;
    for i = 1:N_vgcc
        z = rand(1,1);
        if (state_PQ(i) == 0)
            if(z < alpha_1PQ*dt)
                state_PQ(i) = 1;
            end
        elseif (state_PQ(i) == 1)
            if(z < alpha_2PQ*dt)
                state_PQ(i) = 2;
            elseif (z >= alpha_2PQ*dt && z < (alpha_2PQ + beta_1PQ)*dt)
                state_PQ(i) = 0;                  
            end
        elseif (state_PQ(i) == 2)
            if(z < alpha_3PQ*dt)
                state_PQ(i) = 3;
            elseif (z >= alpha_3PQ*dt && z < (alpha_3PQ + beta_2PQ)*dt)
                state_PQ(i) = 1;                  
            end
        elseif (state_PQ(i) == 3)
            if(z < alpha_4PQ*dt)
                state_PQ(i) = 4;
            elseif (z >= alpha_4PQ*dt && z < (alpha_4PQ + beta_3PQ)*dt)
                state_PQ(i) = 2;                  
            end       
        elseif (state_PQ(i) == 4)
            if(z < beta_4PQ*dt)
                state_PQ(i) = 3;
            end                  
        end 
        
        if(state_PQ(i) == 4)
            OpenPQ = OpenPQ + 1;
        end
        
    end  
          
    PoPQ = OpenPQ/N_vgcc;
    
    IcaPQ = gPQ*PoPQ*(V - Eca);                % uA

    
end 

