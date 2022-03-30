function dydt = ODEs(t, y, Po_ip3r, IcaPQ, AP_condition, cell_condition, ISI, abeta_dose, coupling_condition)

global N_vgcc

%% Solution to ODE's at previous time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cac = y(1);              % Cytosolic Calcium concentration 
Ca_ipr = y(2);           % IP3R microdomain Calcium concentration 
Ct  = y(3);              % Total intracellular Calcium concentration  
V = y(4);                % Membrane potential
hNa = y(5);              % Sodium current inactivating gating variable 
nK = y(6);               % Potassium current activating gating variable 
PLC = y(7);              % PLC Concentration
G = y(8);                % G-proteins
IP3 = y(9);              % IP3 concentration
Glut = y(10);            % Glutamate concentration %% ''' not used ''' %%
Ca_vgcc = y(11);         % VGCC nanodomain calcium concentration 

%% Variables for Allosteric Model

% Vesicle number in: Reserve pool - R_allo, Docked pool - U_allo,
% Vesicle number in: Slow Release Pool (SRP) Refractory pool - RFv_allo, Flow Release Pool (FRP) Refractory pool - RFv_allo 
% Vesicle number in: V(0-1) - SRP, W(0-1) - FRP,

R_allo = y(12); U_allo = y(13); RFv_allo = y(14); RFw_allo = y(15);  
V0 = y(16); V1 = y(17); V2 = y(18); V3 = y(19); V4 = y(20); V5 = y(21);
W0 = y(22); W1 = y(23); W2 = y(24); W3 = y(25); W4 = y(26); W5 = y(27);


%% Variables forDual Sensor Model 

% Vesicle number in: Reserve pool - R_dual, Docked pool - U_dual,
% Vesicle number in: Slow Release Pool (SRP) Refractory pool - RFv_dual, Flow Release Pool (FRP) Refractory pool - RFv_dual 
% Vesicle number in: V(0,0 - 5,2) - SRP, W(0,0 - 5,2) - FRP,


R_dual = y(28) ; U_dual = y(29); RFv_dual = y(30); RFw_dual = y(31);

V00 = y(32); V01 = y(33); V02 = y(34); V10 = y(35); V11 = y(36); V12 = y(37); V20 = y(38); V21 = y(39); V22 = y(40); 
V30 = y(41); V31 = y(42); V32 = y(43); V40 = y(44); V41 = y(45); V42 = y(46); V50 = y(47); V51 = y(48); V52 = y(49);

W00 = y(50); W01 = y(51); W02 = y(52); W10 = y(53); W11 = y(54); W12 = y(55); W20 = y(56); W21 = y(57); W22 = y(58); 
W30 = y(59); W31 = y(60); W32 = y(61); W40 = y(62); W41 = y(63); W42 = y(64); W50 = y(65); W51 = y(66); W52 = y(67);


%% Parameter Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    gamma_1 = 100;                      % cytoplasmic-to-ER-microdomain Volume ratio
    gamma_2 = 10;                       % cytoplasmic-to-ER Volume ratio 
    gamma_3 = 60;                       % cytoplasmic-to-VGCC-nanodomain ratio  
    Vol_er = 3.9*0.1*0.1e-18;           % (um)^3 ER volume %% '''not used'''
   
%% Algorithm for applied stimulation

% Applied/Injected Current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that the amplitude of applied current is what sets the frequency of the Action Potential (AP) in the train 

    if AP_condition=="Single"     
        % ''' Apply stimulation only after 3 ms. This is chosen so that the system relaxes before stimulation is applied ''' 
        if t > 3
           Iapp = 0.0;         % uA/cm^2
        else
           Iapp = 10;          % uA/cm^2    
        end
    elseif AP_condition=="Train"
        % ''' Apply stimulation only for roughly about 400 ms ''' 
        if t < 3 || t > 435
           Iapp = 0.0;         % uA/cm^2
        else
           Iapp = 1.7;          % uA/cm^2    Chosen to get a AP frequency of 20 HZ as desired in simulations 
                                %            and can be adjusted to get other desired frequencies
        end
    elseif AP_condition=="Paired Pulse"
        % ''' Employs specified Inter-Spike Interval (ISI) for paired pulse protocol. '''
        % ''' Default ISI = 40 ms''' 
        if t < 3 || t > 3 + ISI
            if t > 3 + 3 + ISI
                Iapp = 0.0;         % uA/cm^2
            else
                Iapp = 10.0;         % uA/cm^2
            end
        else 
            Iapp = 0.0;         % uA/cm^2
        end
    end


    e:
    
%%  Currents

   [INa, IK, IL, dhNadt, dnKdt] = CurrentsandGatingVariables(V, hNa, nK);
        
%% Calculating Endoplasmic Reticulum Calcium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Cer = gamma_2*(Ct - Cac - Ca_ipr/gamma_1 - Ca_vgcc/gamma_3);
    
%% Fluxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PMCA FLUX
    
    Vpmca = 3.195;   % uM/ms   Maximum capacity of plasma pump
    Kpmca = 0.5;     % uM      Half-maximal activating Cac of plasma pump     
    nP = 2.0;        %         Hill coefficient of plasma pump

    Jpmca = Vpmca*(Cac^nP)/(Kpmca^nP + Cac^nP);           % uM/ms  Flux through plasma pump
    
%% SERCA Flux 

    Vs = 10.0;       % uM/ms   Maximum capacity of SERCA
    Ks = 0.26;       % uM      SERCA half-maximal activating Cac
    ns = 1.75;       %         Hill coefficient of SERCA
    
    Jserca = Vs * (Cac^ns) / (Ks^ns + Cac^ns);            % uM/ms  Uptake of Ca into the ER by SERCA pumps

%% VGCC Flux
    

%   Parameters from Joe Latulippe. et. al 2021
    z = 2;                                      %         valence of Ca ion
    F = 96485.33;                               % C/mole  Faraday's constant, unit: coul/mole
    Area = 3.8489e-09;                          % cm^2    Bouton Membrane area in cm^2
    Vol_tmnal = 1.22e-16;                       % Litre   Presynaptic terminal Volume (Assuming a spherical Bouton shape)
    Kdiff_vgcc = 0.071;                         % /ms     Rate of diffusion from VGCC nanodomain
    cluster_radius = 25e-03;                    % um
    cluster_area = (pi)*(cluster_radius)^2;     % um^2
    active_zone_area = 0.04;                    % um^2
    active_zone_number = 1.3;                   %         Number of active zones
    
    channel_density = N_vgcc / (active_zone_number * active_zone_area);
    
    Ica = channel_density * cluster_area * IcaPQ;  
    
    % Macroscopic current density through an ensemble of open channels where channel_density is crucial:
    Jvgcc = (-Ica/(z*F*Vol_tmnal))*1e-03;         % uM/ms   Calcium influx from the extracellular space to cytosol(with Mitoc 2.03)
    Jdiff_vgcc = Kdiff_vgcc * (Ca_vgcc - Cac);    % uM/ms   Diffusion from VGCC cluster nanodomain to the cytoplasm
    
%% Receptors -- IP3 and RYR

    KIPR = 5;                                      % /ms  IP3R flux coefficient 
    kdiff = 10;                                    % /ms  Ca diffusional flux coefficient
    
    Jipr  = KIPR * Po_ip3r * (Cer - Ca_ipr);       % uM/ms  Flux through the IP3R
    Jdiff = kdiff * (Ca_ipr - Cac);                % uM/ms  Diffusion from ER cluster nanodomain to the cytoplasm
    
 %% Leak Fluxes -- ER and Plasma Membrane

    kleak = 0.0022;                   % |   /ms     |	ER leak flux coefficient
    Jleakin = 0.03115;                % |   uM/ms   |	Plasma membrane leak influx
    Vleakin = 0.2;                    % |   /ms     |	IP3 In Leak flux coefficient

    Jleak = kleak * (Cer - Cac);      % |   uM/ms   |	Baground leak from the ER into the cytoplasm    
    Jin = Jleakin + Vleakin*IP3; 

 %% Flux from VGCC Nano Domain to IP3R Nanodomain
    if coupling_condition == "Higher_Coupling_AD"
        if cell_condition == "WT"
            Kc = 20; 
            Vc = 118;
            K_hat = 5;
        elseif cell_condition == "AD"
            Kc = 10; 
            Vc = 118;
            K_hat = 15;
        end
    elseif coupling_condition == "Higher_Coupling_WT"
        if cell_condition == "WT"
            Kc = 10; 
            Vc = 118;
            K_hat = 15;
        elseif cell_condition == "AD"
            Kc = 20; 
            Vc = 118;
            K_hat = 5;
        end
    elseif coupling_condition == "Same_Coupling"
            Kc = 10; 
            Vc = 118;
            K_hat = 15;
    end
    
    JCa_vgcc_Ca_ipr = Vc * (Ca_vgcc^2 - K_hat*Ca_ipr^2)/(Kc^2 + Ca_vgcc^2);
    
        
%% Differential Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Membrane Potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    C = 1;                                                  % uF/cm^2  membrane capacitance
    dvdt = (0*Iapp + INa + IK + IL)/C;   

% Calcium Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dCacdt  =  (Jdiff_vgcc + Jin + Jdiff + Jleak - Jpmca - Jserca);
    dCa_iprdt  =  gamma_1*(Jipr - Jdiff) + JCa_vgcc_Ca_ipr;
    dCa_vgccdt  =  gamma_3*(Jvgcc - Jdiff_vgcc - JCa_vgcc_Ca_ipr/gamma_1);
    dCtdt  = (Jin + Jvgcc - Jpmca);

%% IP3 Production %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          
    
    % Parameters from Joe Latulippe. et. al 2021
    % Contains Amyloid beta effect.
    
    G_tot = 1;                          %           Scaled total number of G-protein
    PLC_tot = 1;                        %           Scaled total number of PLC
    rbeta = 0.001;                      %           A? decay rate
    Kq = 0.0086;                        % ug/mL     PLC dissociation constant
    
    if cell_condition == "WT"    
        
        abeta = 0;                          %           Amyloid beta concentration
        Vo = 0.15;                          % uM        Intrinsic PLC-mediated IP3 production
        Vq = 7.82;                          % uM        Control parameter for influence of A? on IP3
        Kip3k = 0.6;                        % uM        Half-activation for 3-kinase
        Kplc = 0.01;                        % uM        PLC sensitivity to Ca2+
        k_3k = 1.5e-03;                     % /ms       IP3 Phosphorylation rate 
        k_5p = 0.01e-03;                    % /ms       IP3 DePhosphorylation rate 
        Kf_plc = 0.35e-03;                  % /ms       PLC-protein activation rate
        Kb_plc = 2.2e-02;                   % /ms       PLC-protein deactivation rate
        Kf_gp = 0.33e-03;                   % /ms       G-protein activation rate
        Kb_gp = 2.17e-03;                   % /ms       G-protein deactivation rate
        delta_G = 0.01;                     %           G-protein intrinsic backgroung activity
        Vr = 7.4;                           %           Maximal G-protein activation
        Kr = 4467;                          % ug/mL     Amyloid beta concentration producing half-activation 
        
    elseif cell_condition == "AD"
        
        abeta = abeta_dose;                %            Amyloid beta concentration       
        Vo = 0.19;                         % uM         Intrinsic PLC-mediated IP3 production
        Vq = 980;                          % uM         Control parameter for influence of A? on IP3
        Kip3k = 1.6;                       % uM         Half-activation for 3-kinase
        Kplc = 0.016;                      % uM         PLC sensitivity to Ca2+
        k_3k = 0.7e-03;                        % /ms    IP3 Phosphorylation rate 
        k_5p = 0.005e-03;                      % /ms    IP3 DePhosphorylation rate 
        Kf_plc = 0.75e-03;                     % /ms    PLC-protein activation rate
        Kb_plc = 2.0e-01;                      % /ms    PLC-protein deactivation rate
        Kf_gp = 0.047e-03;                     % /ms    G-protein activation rate
        Kb_gp = 4.7e-03;                       % /ms    G-protein deactivation rate
        delta_G = 0.012;                   %            G-protein intrinsic backgroung activity
        Vr = 10;                           %            Maximal G-protein activation
        Kr = 2000;                         % ug/mL      Amyloid beta concentration producing half-activation
        
    end
    
    dGlutdt = 0;    % Not used! 
    
    t_ip3 = 1/(k_3k + k_5p);            % ms            represents the characteristic time of IP3 turnover
    Nu = k_3k/(k_3k + k_5p);            %               affects negative feedback from Ca2+ to IP3 degredation
    q = heaviside(t - 2) * abeta * exp(-rbeta * (t - 2) * heaviside(t - 2));
    rho_G = Vr * q/(Kr + q);                             % Amyloid beta dependent G-protein activation
    Vplc = Vo + (Vq * q^2/(Kq^2 + q^2));                 % Amyloid beta dependent maximal rate of PLC mediated IP3 production                                 
    Jplc = Vplc * PLC * (Cac^2 / (Kplc^2 + Cac^2));
    Jdeg = (Nu * (Cac^2 / (Kip3k^2 + Cac^2)) + (1 - Nu)) * IP3;
   
    dPLCdt = Kf_plc * G * (PLC_tot - PLC) - Kb_plc * PLC;
    dGdt = Kf_gp * (rho_G + delta_G)*(G_tot - G) - Kb_gp*G;
    dIP3dt = (Jplc - Jdeg)/t_ip3;

    
%% Transmitter Release Rate Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Allosteric Model

  
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

    

%% Return Solutions
  
% Deterministic Solutions

dydt = [dCacdt; dCa_iprdt; dCtdt; dvdt; dhNadt; dnKdt;...
        dPLCdt; dGdt; dIP3dt; dGlutdt;...
        dCa_vgccdt;...
        dR_allodt; dU_allodt; dRFv_allodt; dRFw_allodt; dV0dt; dV1dt; dV2dt; dV3dt; dV4dt; dV5dt; dW0dt; dW1dt; dW2dt;...
        dW3dt; dW4dt; dW5dt;...
        dR_dualdt; dU_dualdt; dRFv_dualdt; dRFw_dualdt;...
        dV00; dV01; dV02; dV10; dV11; dV12; dV20; dV21; dV22; dV30; dV31; dV32; dV40; dV41; dV42; dV50; dV51; dV52;...
        dW00; dW01; dW02; dW10; dW11; dW12; dW20; dW21; dW22; dW30; dW31; dW32; dW40; dW41; dW42; dW50; dW51; dW52;...
        ];
