function dydt = ODEs(t, y, Po_ip3r, IcaPQ, AP_condition, cell_condition, ISI, abeta_dose, coupling_condition, pore_open_times, calcium_source)

global N_pores N_vgcc

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
Cabeta = y(10);          % Calcium concentration in the sub-plasmalemmal compartment of the Amyloid beta pore
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
    gamma_4 = 100;                      % cytoplasmic-to-Abeta-subplasmalemmal compartment Volume ratio  
   
%% Parameter Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    gamma_1 = 100;                      % cytoplasmic-to-ER-microdomain Volume ratio
    gamma_2 = 10;                       % cytoplasmic-to-ER Volume ratio 
    gamma_3 = 60;                       % cytoplasmic-to-VGCC-nanodomain ratio  
    Vol_er = 3.9*0.1*0.1e-18;           % (um)^3 ER volume %% '''not used'''
   
%% Algorithm for applied stimulation

% Applied/Injected Current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that the amplitude of applied current is what sets the frequency of the Action Potential (AP) in the train 

    [Iapp]= CurrentInjection(t, AP_condition, ISI);
    
%%  Currents

   [INa, IK, IL, dhNadt, dnKdt] = CurrentsandGatingVariables(Cac, V, hNa, nK);
        
%% Calculating Endoplasmic Reticulum Calcium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Cer = gamma_2*(Ct - Cac - Ca_ipr/gamma_1 - Ca_vgcc/gamma_3);


%% Fluxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Ica_density, Iabeta_density, Jpmca, Jserca, Jvgcc, Jdiff_vgcc, Jipr, Jdiff,...
    Jleak, Jin, JCa_vgcc_Ca_ipr, Jabeta, Jdiff_abeta,...
    Jdiff_abeta_cyt] = Fluxes(V, Cac, Ca_vgcc, Cabeta, Cer, Ca_ipr, IP3, ...
                             IcaPQ, Po_ip3r, coupling_condition, cell_condition,...
                              t, pore_open_times, calcium_source);

        
%% Differential Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Membrane Potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    C = 1;                                                 % uF/cm^2  membrane capacitance
    dvdt = (Iapp + INa + IK + IL + Ica_density + Iabeta_density)/C;  

%% Calcium Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
    dCacdt  =  (Jdiff_vgcc + Jin + Jdiff + Jleak - Jpmca - Jserca);
    dCa_iprdt  = gamma_1*(Jipr - Jdiff) + JCa_vgcc_Ca_ipr;
    dCa_vgccdt = gamma_3*(Jvgcc + Jabeta - Jdiff_vgcc - JCa_vgcc_Ca_ipr/gamma_1);
    dCabetadt = gamma_4*(Jabeta - Jdiff_abeta/gamma_3 - Jdiff_abeta_cyt);
    dCtdt  = (Jin + Jvgcc + Jabeta - Jpmca);

%% IP3 Production %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          
    
    [Jplc, Jdeg, t_ip3, dPLCdt, dGdt, dGlutdt] = IP3Production(t, Cac, IP3, PLC, G, abeta_dose, cell_condition);
    
    dIP3dt = (Jplc - Jdeg)/t_ip3;

    
%% Transmitter Release Rate Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Allosteric Model

    [dR_allodt, dU_allodt, dRFv_allodt, dRFw_allodt, dV0dt,...
          dV1dt, dV2dt, dV3dt, dV4dt, dV5dt, dW0dt, dW1dt, dW2dt,...
          dW3dt, dW4dt, dW5dt] = AllostericModel(Cac, Ca_vgcc, R_allo,...
          U_allo, RFv_allo, RFw_allo, V0, V1, V2, V3, V4, V5, W0, W1, W2,...
          W3, W4, W5);

    
%% Dual Sensor Model

 [...
        dR_dualdt, dU_dualdt, dRFv_dualdt, dRFw_dualdt, dV00, dV01, dV02, dV10, ...
        dV11, dV12, dV20, dV21, dV22, dV30, dV31, dV32, dV40, dV41, dV42, dV50,...
        dV51, dV52, dW00, dW01, dW02, dW10, dW11, dW12, dW20, dW21, dW22, dW30,...
        dW31, dW32, dW40, dW41, dW42, dW50, dW51, dW52...
         ] = DualSensorModel(Cac, Ca_vgcc, R_dual, U_dual, RFv_dual, RFw_dual, ...
                            V00, V01, V02, V10, V11, V12, V20, V21, V22, V30, ...
                            V31, V32, V40, V41, V42, V50, V51, V52, W00, W01,...
                            W02, W10, W11, W12, W20, W21, W22, W30, W31, W32, ...
                            W40, W41, W42, W50, W51, W52);

%% Return Solutions
  
% Deterministic Solutions

dydt = [dCacdt; dCa_iprdt; dCtdt;... 
        dvdt; dhNadt; dnKdt;...
        dPLCdt; dGdt; dIP3dt; dCabetadt;...
        dCa_vgccdt;...
        dR_allodt; dU_allodt; dRFv_allodt; dRFw_allodt; dV0dt; dV1dt; dV2dt; dV3dt; dV4dt; dV5dt; dW0dt; dW1dt; dW2dt;...
        dW3dt; dW4dt; dW5dt;...
        dR_dualdt; dU_dualdt; dRFv_dualdt; dRFw_dualdt;...
        dV00; dV01; dV02; dV10; dV11; dV12; dV20; dV21; dV22; dV30; dV31; dV32; dV40; dV41; dV42; dV50; dV51; dV52;...
        dW00; dW01; dW02; dW10; dW11; dW12; dW20; dW21; dW22; dW30; dW31; dW32; dW40; dW41; dW42; dW50; dW51; dW52;...
        ];

