function    [Ica_density, Iabeta_density, Jpmca, Jserca, Jvgcc, Jdiff_vgcc, Jipr, Jdiff,...
            Jleak, Jin, JCa_vgcc_Ca_ipr, Jabeta, Jdiff_abeta,...
            Jdiff_abeta_cyt] = Fluxes(V, Cac, Ca_vgcc, Cabeta, Cer, Ca_ipr, IP3, ...
                                      IcaPQ, Po_ip3r, coupling_condition, cell_condition,...
                                      time, pore_open_times, calcium_source)

global N_vgcc N_pores SCL

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
    
%   Parameters (some modified) from Schikorski. et.
%   al. 1997, Holderith et. al. 2012, and Ermolyuk. et. al 2013, 

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
    Ica_density = channel_density * IcaPQ;

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

    if cell_condition == "WT"
        Kc = 20; 
        Vc = 118;
        K_hat = 5;
    elseif cell_condition == "AD"
        Kc = 10; 
        Vc = 118;
        K_hat = 15;
    end
    
    JCa_vgcc_Ca_ipr = Vc * (Ca_vgcc^2 - K_hat*Ca_ipr^2)/(Kc^2 + Ca_vgcc^2);
    
%% Flux through Amyloid Beta Pores

    % Parameters from Prista von Bonhorst, F.; Gall, D.; Dupont, G. 
    % Impact of β-Amyloids Induced Disruption of Ca2+ Homeostasis in a Simple Model of Neuronal Activity.
    % Cells 2022, 11, 615. https://doi.org/10.3390/cells11040615
    
    Vabeta = 2;                                                       % uM/ms     Maximal rate of Ca2+ entry through the pores ()    
    buffered_diff_const = 0.015;                                      % (um^2)/ms Buffered diffusion coefficient (Allbritton, et. al 1992 Science)                             
    vgcc_diff_distance = 0.035;                                       % um        Diffusional distance to VGGC nanodomain
    cyto_diff_distance = 0.35;                                        % um        Diffusional distance to bulk cytosol
    Kdiff_abeta_vgcc = buffered_diff_const/(vgcc_diff_distance^2);    % /ms       Abeta-VGCC nanodomain diffusion coefficient. 
    Kdiff_abeta_cyto = buffered_diff_const/(cyto_diff_distance^2);    % /ms       Abeta-Cytosol nanodomain diffusion coefficient. 

    % Parameters Syed. Et al. Modeling the kinetics of amyloid beta pores
    % and long-term evolution of their Ca2+ toxicity 2022
    % This version posted May 3, 2022.; https://doi.org/10.1101/2022.05.02.490365


    rpore = 2.6e-03;                                        % um      Radius of pore (Ref. At least 1.3nm from Lin H. 2001, Jang H 2007., Sepulveda. 2010)
    pore_current = SCL*0.5e-06;                                 % uA      Current per permeability level Demuro et al 2011 J Cell Biol
    pore_domain_number = 1;                                 %         Number of independent clusters. 
    pore_domain_radius = 0.08;                              % um      Raidus of sub-plasmalemmal region/domain of calcium around pore.
    pore_domain_area = pi*pore_domain_radius^2;             % um^2    Area of sub-plasmalemmal region/domain of calcium around pore. Syed et. al version posted April 30, 2022. ; https://doi.org/10.1101/2022.04.29.490101
    pore_domain_volume = (2/3)*pi*pore_domain_radius^3;     % um^3    Volume of sub-plasmalemmal region/domain of calcium around pore.
    pore_radius = rpore;                                    % um
    Vol_pore = (1e-15)*(2/3)*pi*(rpore^3);                  % Litre   Volume of hemispher3 over the pore (Assuming a spherical pore/channel)
    pore_cluster_area = (pi)*(pore_radius)^2;               % um^2    Area of cluster of pores. 

    pore_density = N_pores / (active_zone_number * active_zone_area);
    
    if any(pore_open_times==time)
        Iabeta = pore_density * pore_cluster_area * pore_current;
        Iabeta_density = pore_density * pore_current;

        %fabeta = Vabeta*(1/(1+ exp((V + 30)/23)));                           % uM/ms   Calcium influx from the extracellular space to cytosol (Prista von Bonhorst, F. et. al. Cells 2022)
        %Iabeta_density = fabeta*2*F*(pore_domain_volume*1e-15)/pore_density;
        
    else
        Iabeta = 0;
        fabeta = 0;
        Iabeta_density = 0;
    end


    if (calcium_source == "ip3r_nc_and_abeta" || calcium_source == "ip3r_hc_and_abeta")
        Jabeta = (Iabeta/(z*F*Vol_tmnal))*1e-03;            % uM/ms   Calcium influx from the extracellular space to cytosol
        %Jabeta = fabeta; 
        Jdiff_abeta = Kdiff_abeta_vgcc * (Cabeta - Ca_vgcc);   % uM/ms   Diffusion from Abeta pore sub-plasmalemmal compartment
                                                               %         to VGCC nanodomain
        Jdiff_abeta_cyt = Kdiff_abeta_cyto * (Cabeta - Cac);   % uM/ms   Diffusion from Abeta pore sub-plasmalemmal compartment
                                                               %         to bulk cytosol         
    else
        Jabeta = 0;
        Jdiff_abeta = 0;
        Jdiff_abeta_cyt = 0;
    end

    

  