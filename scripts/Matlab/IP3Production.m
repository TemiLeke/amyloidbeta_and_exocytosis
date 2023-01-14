function  [Jplc, Jdeg, t_ip3, dPLCdt, dGdt, dGlutdt] = IP3Production(t, Cac, IP3, PLC, G, abeta_dose, cell_condition)


%% IP3 Production %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parameters (some modified) from Joe Latulippe. et. al 2021
    % Contains Amyloid beta effect.

    rbeta = 0.001;                      %           A? decay rate
    Kq = 0.0086;                        % ug/mL     PLC dissociation constant
    G_tot = 1;                          %           Scaled total number of G-protein
    PLC_tot = 1;                        %           Scaled total number of PLC
    
    if cell_condition == "WT"    
        
        abeta = 0;                          %           Amyloid beta concentration
        Vo = 0.15;                          % uM        Intrinsic PLC-mediated IP3 production
        Vq = 7.82;                          % uM        Control parameter for influence of A? on IP3
        Kip3k = 0.6;                        % uM        Half-activation for 3-kinase
        Kplc = 0.01;                        % uM        PLC sensitivity to Ca2+
        k_3k = 1.5e-03;                     % /ms       IP3 Phosphorylation rate 
        k_5p = 0.01e-03;                    % /ms       IP3 DePhosphorylation rate 
        Vr = 7.4;                           %           Maximal G-protein activation
        Kr = 4467;                          % ug/mL     Amyloid beta concentration producing half-activation 
        Kf_plc = 0.35e-03;                  % /ms       PLC-protein activation rate
        Kb_plc = 2.2e-02;                   % /ms       PLC-protein deactivation rate
        Kf_gp = 0.33e-03;                   % /ms       G-protein activation rate
        Kb_gp = 2.17e-03;                   % /ms       G-protein deactivation rate
        delta_G = 0.01;                     %           G-protein intrinsic backgroung activity  
        
        
    elseif cell_condition == "AD"
        
        abeta = abeta_dose;                 %            Amyloid beta concentration       
        Vo = 0.19;                          % uM         Intrinsic PLC-mediated IP3 production
        Vq = 980;                           % uM         Control parameter for influence of A? on IP3
        Kip3k = 1.6;                        % uM         Half-activation for 3-kinase
        Kplc = 0.016;                       % uM         PLC sensitivity to Ca2+
        k_3k = 0.7e-03;                     % /ms        IP3 Phosphorylation rate 
        k_5p = 0.005e-03;                   % /ms        IP3 DePhosphorylation rate 
        Vr = 10;                            %            Maximal G-protein activation
        Kr = 2000;                          % ug/mL      Amyloid beta concentration producing half-activation
        Kf_plc = 0.75e-03;                  % /ms       PLC-protein activation rate
        Kb_plc = 2.0e-01;                   % /ms       PLC-protein deactivation rate
        Kf_gp = 0.047e-03;                  % /ms       G-protein activation rate
        Kb_gp = 4.7e-03;                    % /ms       G-protein deactivation rate
        delta_G = 0.012;                    %           G-protein intrinsic backgroung activity
        
    end
    
    t_ip3 = 1/(k_3k + k_5p);            % ms            represents the characteristic time of IP3 turnover
    Nu = k_3k/(k_3k + k_5p);            %               affects negative feedback from Ca2+ to IP3 degredation
    q = heaviside(t - 2) * abeta * exp(-rbeta * (t - 2) * heaviside(t - 2));
    rho_G = Vr * q/(Kr + q);                             % Amyloid beta dependent G-protein activation
    Vplc = Vo + (Vq * q^2/(Kq^2 + q^2));                 % Amyloid beta dependent maximal rate of PLC mediated IP3 production                                 
    Jplc = Vplc * PLC * (Cac^2 / (Kplc^2 + Cac^2));
    Jdeg = (Nu * (Cac^2 / (Kip3k^2 + Cac^2)) + (1 - Nu)) * IP3;
    
    dGlutdt = 0;    % Not Used!
    dPLCdt = Kf_plc * G * (PLC_tot - PLC) - Kb_plc * PLC;
    dGdt = Kf_gp * (rho_G + delta_G)*(G_tot - G) - Kb_gp*G;
        
