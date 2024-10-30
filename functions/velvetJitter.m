% Function for generating jittered 2-channel velvet noise

function [seq,k,sign] = velvetJitter(len, Td, decay, jitter, jitterAlg, varargin)
    %Td = fs/density;                               % average spacing between impulses = grid size
    M = round(len/Td);                              % Total number of impulses in the sequence
    if numel(jitter)~=M && ~isempty(jitter)
        error("Length of jitter vector does not match the number of pulses!")
    end
    if sum(jitter > Td) ~= 0
        warning(['Jitter larger than grid size for ' num2str(sum(jitter > Td)), '/', num2str(M), ' pulses.'])
    end
    r1 = rand(M,1);                                 % random number sequence
    r2 = rand(M,1);                                 % random number sequence
    r3 = rand(M,1);                                 % random number sequence
    
    alpha = (-log(10.^(-decay/20)))/len;            % decay constant
    s = 2 * round(r1) - 1;                          % signs
    k = zeros(M,2);
    if isempty(jitter)
        k(:,1) = ceil((0:M-1)'.*Td + r2.*(Td-1)); 
        k(:,2) = ceil((0:M-1)'.*Td + r3.*(Td-1));
    else
    switch jitterAlg
        case 'Additive'
            k(:,1) = ceil((0:M-1)'.*Td + r2.*(Td-1));       % pulse locations left channel
            k(:,2) = max(1, k(:,1) + jitter);               % pulse locations right channel
        case 'Replacing'
            k(:,1) = ceil((0:M-1)'.*Td + r2.*(Td-1));       % pulse locations left channel
            k(:,2) = max(1, k(:,1) + jitter);               % pulse locations right channel
        case 'GridLock'
            for m = 1:M
                if (jitter(m)>=Td)                          % out of bounds positive jitter                           
                    k(m,1) = ceil((m-1).*Td+eps);      
                    k(m,2) = ceil((m-1).*Td + (Td-1)); 
                elseif (jitter(m)<=-Td)                     % out of bounds negative jitter
                    k(m,1) = ceil((m-1).*Td + (Td-1));      
                    k(m,2) = ceil((m-1).*Td+eps); 
                elseif (jitter(m)<Td)                       % within bounds positive jitter 
                    k(m,1) = ceil((m-1)*Td + r2(m)*(Td-1-jitter(m)));
                    k(m,2) = max(1,k(m,1) + jitter(m));
                elseif (0>jitter(m)&&jitter(m)>-Td)         % within bounds negative jitter
                    k(m,1) = ceil((m-1)*Td + Td+jitter(m)+r2(m)*(abs(jitter)-1));
                    k(m,2) = max(1,k(m,1) + jitter(m));
                    
                end
            end
    end
    end
    %k = min(len,k); % limit to the desired length
    
    seq = zeros(max(len,max(k, [], 'all')),2);
    
    seq(k(:,1),1) = s .* exp(-alpha*k(:,1));
    sign(:,1) = s;
    if numel(varargin)~=0 && varargin == "independentSign"
        s = 2 * round(rand(M,1)) - 1;   % update signs for the other channel   
    end
    switch jitterAlg
        case 'Additive'
            for m = 1:M
                seq(k(m,2),2) = seq(k(m,2),2) + s(m) .* exp(-alpha*k(m,2));
            end
        case 'Replacing'
            seq(k(:,2),2) = s .* exp(-alpha*k(:,2));
        case 'GridLock'
            seq(k(:,2),2) = s .* exp(-alpha*k(:,2));
    end
    seq = seq(1:len,:); % enforce length
    sign(:,2) = s;
end