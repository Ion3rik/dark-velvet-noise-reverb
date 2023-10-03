% Function for generating time-varying Velvet Pulse Locations
% 9.2.2023
% Jon Fagerström

% Arguments:
%   <fs>:   Sampling frequency
%   <len>:  Length of the generated DVN sequence
%   <rho>:  Density pulses/s [start, end]

% Returns:
%       <k> pulse locations with shape (M,1)
%       <s> pulse signs with shape (M,1)

function [k,s, seq] = vnTimeVar(fs, rhoRange, len)
    % TIME VARYING PARAMETERS %
    %TdRange = fs ./ rhoRange; % average spacing between pulses = grid size
    
    rho = getTimeVaryingParameterValue(rhoRange, len, 'Linear');
    Td = fs ./ rho;
   
    % PULSE LOCATIONS %
    pulseIdx = 1;
    m = 2;
    k(1) = 1;
    Tdcum = Td(1);
    while pulseIdx <= len
        r2 = rand;
        n = k(m-1); % previous pulse location
        pulseIdx = round(Tdcum + r2*(Td(n)-1));
        Tdcum = Tdcum + Td(n);
        k(m) = pulseIdx;
        m = m+1;
    end
    k = k(1:end-1);
    M = length(k); % get the total amount of pulses
    
    % PULSE SIGNS %
    r1 = rand(1,M); % random number sequence for the sign
    s(1,1:M) = 2*round(r1(1,1:M))-1; 

    % GENERATE THE SEQUENCE
    seq = zeros(len,1);
    for i = 1:length(s)                    % sequence
        s(i) = s(i) * sqrt(Td(k(i)));             % compensate the varibale density
        seq(k(i)) = s(i);
    end
    seq = seq(1:len);
end

function modEnv = getTimeVaryingParameterValue(range, len, type)
    switch type
        case 'Linear'
            modEnv = linspace(range(1), range(2), len);
        case 'Log'
            modEnv = log10(linspace(1,10,len)) * (range(2)-range(1)) + range(1);
%             modEnv = exp(linspace(log(range(1)), log(range(2)),len));
            
    end
end

