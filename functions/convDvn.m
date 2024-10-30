% Function for implementing the DVN convolution
% Jon Fagerström
% 23.9.2022
% Arguments:
%           <x>: mono input audio, with size (numSamples,1)
%           <k>: pulse locations, with size (numPulses,1)
%           <g>: pulse gains, with size (numPulses,1)
%           <B>: pulse filter coeffs, with size (feedforward order, numFilters, numStages)
%           <A>: pulse filter coeffs, with size (feedback order, numFilters, numStages)
%           <gAll>: Allpass filter gains
%           <dAll>: Allpass filter delays
%           <filtList>: list of which pulses connects to which filter, with size (numPulses,1)
% Outputs:
%           <y>: DVN convolution result    
function [y, filterOut, whiteOut] = convDvn(x,k,g,filtList,B, A)
    numFilters = size(A,2); % number of filters
    vL = max(k);            % maximum delay of the pulses
    numStages = size(B,3);   % number of stages in the pulse filters
    
    % Split processing to each filter branch
    filterOut = zeros(length(x) + vL - 1, numFilters);
    whiteOut = zeros(length(x) + vL - 1, numFilters); 
    for f = 1:numFilters
        kf = k(filtList == f);                                      % pick each pulse location that routes to the current filter branch
        gf = g(filtList == f);                                      % pick each pulse gain that routes to the current filter branch
        velvetOut = convSparse(x,kf,gf,vL);                         % convolve input with the velvet taps of this filter branch
        for s = 1:numStages
            filterOut(:,f) = filter(B(:,f,s),A(:,f,s), velvetOut);      % convolve velvetOut with the pulse filters, and add to output
            whiteOut(:,f) = velvetOut;
        end
    end
    y = sum(filterOut,2);
end

