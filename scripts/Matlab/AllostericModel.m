function [dR_allodt, dU_allodt, dRFv_allodt, dRFw_allodt, dV0dt,...
          dV1dt, dV2dt, dV3dt, dV4dt, dV5dt, dW0dt, dW1dt, dW2dt,...
          dW3dt, dW4dt, dW5dt] = AllostericModel(Cac, Ca_vgcc, R_allo,...
          U_allo, RFv_allo, RFw_allo, V0, V1, V2, V3, V4, V5, W0, W1, W2,...
          W3, W4, W5)

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
    