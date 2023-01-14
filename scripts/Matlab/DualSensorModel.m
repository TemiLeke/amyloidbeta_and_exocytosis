function [...
        dR_dualdt, dU_dualdt, dRFv_dualdt, dRFw_dualdt, dV00, dV01, dV02, dV10, ...
        dV11, dV12, dV20, dV21, dV22, dV30, dV31, dV32, dV40, dV41, dV42, dV50,...
        dV51, dV52, dW00, dW01, dW02, dW10, dW11, dW12, dW20, dW21, dW22, dW30,...
        dW31, dW32, dW40, dW41, dW42, dW50, dW51, dW52...
         ] = DualSensorModel(Cac, Ca_vgcc, R_dual, U_dual, RFv_dual, RFw_dual, ...
                            V00, V01, V02, V10, V11, V12, V20, V21, V22, V30, ...
                            V31, V32, V40, V41, V42, V50, V51, V52, W00, W01,...
                            W02, W10, W11, W12, W20, W21, W22, W30, W31, W32, ...
                            W40, W41, W42, W50, W51, W52)

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

    