function [INa, IK, ICl, dhNadt, dnKdt] = CurrentsandGatingVariables(Cac, V, hNa, nK)


%% Nernst potential for Na, K, and Cl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ENa = 50;    %  mV                           
    EK = -100;   %  mV                          
    ECl = -70;   %  mV    
    
%%  Currents

    gNa = 120;           % mS/cm^2  Sodium conductane
    gNaLeak = 0.0175;    % mS/cm^2  Sodium Leak conductane
    gK = 36;             % mS/cm^2  Potassium Leak conductane
    gKLeak = 0.05;       % mS/cm^2  Potassium Leak conductane
    gClLeak = 0.05;      % mS/cm^2  Chlorine Leak conductane
    phi = 5.0;
    gAHP = 0.01;  
                 
    % Gating Variables For Current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    alphan = 0.01*(V + 34)/(1 - exp(-0.1*(V + 34)));        % /ms    
    betan = 0.125*exp(-(V + 44)/80);                        % /ms
    alphah = 0.07*exp(-(V + 44)/20);                        % /ms
    betah = 1.0/(1 + exp(- 0.1*(V + 14)));                  % /ms
    alpham = 0.1*(V + 30)/(1 - exp(-0.1*(V + 30)));         % /ms
    betam = 4*exp(-(V + 55)/18);                            % /ms
    m1 = alpham*V/(alpham*V + betam*V);
                                     
    INa = -gNa*(m1^3)*hNa*(V - ENa)- gNaLeak*(V - ENa);                   % uA/cm^2
    IK  = -(gK*nK^4 + (gAHP*Cac)/(1 + Cac))*(V - EK) - gKLeak*(V - EK);  % uA/cm^2
    ICl = -gClLeak*(V - ECl);                                            % uA/cm^2
 
    
    dnKdt = phi*(alphan*(1 - nK) - betan*nK);
    dhNadt = phi*(alphah*(1 - hNa) - betah*hNa);