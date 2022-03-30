function [INa, IK, IL, dhNadt, dnKdt] = CurrentsandGatingVariables(V, hNa, nK)

%{
    This model is a previously published model of hippocampal pyramidal
    neurons by Kopell et al., 2010 [https://doi.org/10.1007/978-1-4419-0996-1_15]

    In this model, The Na+ and K+ currents consist of active currents 
    corresponding to voltage gated Na+ channels (VGSC) and  delayed
    rectifier K+ channels , passive leak currents

%}


%% Transient NaT current; from Kopell et al., 2010

    V_mNa_th_alpha = 54;     %   mV      | Half activating potential for forward rate(i.e midpoint of the activation curve of NaT 
    k_mNa_alpha = 4.0;       %   mV      | Activation slope for forward rate
    V_mNa_th_beta = 27;      %   mV      | Half activating potential for backward rate(i.e midpoint of the activation curve of NaT 
    k_mNa_beta = 5.0;        %   mV      | Activation slope for backward rate       
   
%%%%%%%%%%%%%%%%%%%% Dynamics of activating gating variable %%%%%%%%%%%%%%%%%%%%%%%%%

    alpha_mNa = 0.32*(V + V_mNa_th_alpha)/(1 - exp(-(V + V_mNa_th_alpha)/k_mNa_alpha));  %   activation variable forward rate
    beta_mNa = -0.28*(V + V_mNa_th_beta)/(1 - exp(-(V + V_mNa_th_beta)/k_mNa_beta));  %   activation variable backward rate
    tau_mNa = 1/(alpha_mNa + beta_mNa);  %   activation variable time constant
    
    mNa_inf = alpha_mNa/(alpha_mNa + beta_mNa);
    
%%%%%%%%%%%%%%%%%%%% Dynamics of inactivating gating variable %%%%%%%%%%%%%%%%%%%%%%%%%    
    

    V_hNa_th_alpha = 50;     %   mV      | Half inactivating potential for forward rate(i.e midpoint of the activation curve of NaT 
    k_hNa_alpha = 18;        %   mV      | inactivation slope for forward rate
    V_hNa_th_beta = 27;      %   mV      | Half inactivating potential for backward rate(i.e midpoint of the activation curve of NaT 
    k_hNa_beta = 5.0;        %   mV      | inactivation slope for backward rate     
       
    
    alpha_hNa = 0.128*exp(-(V + V_hNa_th_alpha)/k_hNa_alpha);  %   inactivation variable forward rate
    beta_hNa =  4.0/(1 + exp(-(V + V_hNa_th_beta)/k_hNa_beta));  %   inactivation variable backward rate
    tau_hNa = 1/(alpha_hNa + beta_hNa);  %   inactivation variable time constant
    
    hNa_inf = alpha_hNa/(alpha_hNa + beta_hNa);
    
    dhNadt = (hNa_inf - hNa)/tau_hNa; 
    
%%%%%%%%%%%%%%%%%%%% Transient NaT current %%%%%%%%%%%%%%%%%%%%%%%%%    

    ENa = 50;                                           %   mV          |   Nernst potential for Na
    gNa = 100;                                          %   mS/cm^2     |   maximal Sodium conductane    
    INa = -gNa*(mNa_inf^3)*hNa*(V - ENa);               %   uA/cm^2     |   Transient sodium cuurent 
    
    
%% The Delayed rectifier K+ current; from Kopell et al., 2010

    V_nK_th_alpha = 52;     %   mV      | Half inactivating potential for forward rate(i.e midpoint of the activation curve of NaT 
    k_nK_alpha = 5.0;       %   mV      | inactivation slope for forward rate
    V_nK_th_beta = 57;      %   mV      | Half inactivating potential for backward rate(i.e midpoint of the activation curve of NaT 
    k_nK_beta = 40.0;       %   mV      | inactivation slope for backward rate     
    

   
    alpha_nK = 0.032*(V + V_nK_th_alpha)/(1 - exp(-(V + V_nK_th_alpha)/k_nK_alpha));  %   inactivation variable forward rate
    beta_nK = 0.5*exp(-(V + V_nK_th_beta)/k_nK_beta);  %   inactivation variable backward rate
    tau_nK = 1/(alpha_nK + beta_nK);  %   inactivation variable time constant
    
    nK_inf = alpha_nK/(alpha_nK + beta_nK);
    
    dnKdt = (nK_inf - nK)/tau_nK;
    
%%%%%%%%%%%%%%%%%%%% Delayed rectifier potassium current %%%%%%%%%%%%%%%%%%%%%%%%%

    EK = -100;                                          %   mV          |   Nernst potential for K+
    gK = 80;                                            %   mS/cm^2     |   maximal Potassium conductance    
    IK = -gK*(nK^4)*(V - EK);                           %   uA/cm^2     |   Deelayed rectifier potassium cuurent 

    
%% The Leak current; from Kopell et al., 2010

    EL = -67;                                           %   mV          |   Nernst potential for Leak current
    gL = 0.1;                                           %   mS/cm^2     |   M-type Potassium conductance    
    IL = -gL*(V - EL);                                  %   uA/cm^2     |   M-type Potassium cuurent 
    
    
%}
       