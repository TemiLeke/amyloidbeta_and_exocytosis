function [Iapp]= CurrentInjection(t, AP_condition, ISI)

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

    