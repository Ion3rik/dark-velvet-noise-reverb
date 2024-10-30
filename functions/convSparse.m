% Function that implements the sparse time domain convolution with a velvet seqeunce
% Author: Jon Fagerström
% Date: 22.9.2022

function y = convSparse(x,k,g,vL)
    xL = length(x);             % length of input signal
    yL = xL + vL - 1;           % length of output vector
    y = zeros(yL,1);            % init output vector
    
    %% SAMPLE BY SAMPLE (Slow!)
    %x = [x; zeros(vL-1,1)];     % zero pad input vector
%     for n = 1:yL
%         pulseMask = (n - k) >= 0;                               % find pulses that contribute to the current output sample
%         y(n) = sum(x(n - k(pulseMask)+1) .* g(pulseMask));      % compute convolution output sample 
%     end
if max(k)>vL
    warning("k-index over the sequence length")
end
    %% PULSE BY PULSE (Still slower than Matlab's conv)
    M = length(k);
    for m = 1:M
        nStart = k(m);
        nStop = nStart + xL - 1;
        y(nStart:nStop) = y(nStart:nStop) + x * g(m);
    end
end