% Function for extracting late reverb from an RIR
% Jon Fagerström
% 28.9.2022
% Arguments:
%           <rir>: room impulse response, with size (numSamples,1)
%           <mixingTime>: time point where to cut the early part in samples
%           <dynamicRange>: dynamic range of the extracted late reverb (cut silent tail out)
% Outputs:
%           <rirLate>: late part of the RIR

function [rirLate, rirEarly] = extractLateReverb(rir, mixingTime, dynamicRange, onset)
    
    if (onset)
        [~,onset] = max(abs(rir));
        rir = rir(onset:end,:);
    end
    %startIdx = find(rir.^2 > 0.001,1);              % find where the signal energy is more -60 dB
    %rirLate = rir(startIdx:end);                    % remove silent begining
    rirEarly = rir(1:mixingTime-1);
    cutPoint = numel(rir);
    if (~isempty(dynamicRange))
        rirEdc = edc(rir);                            % compute the energy decay
        cutPoint = find(rirEdc <= -dynamicRange,1);   % find where to cut the tail based on the dynamic range
    end
    rirLate = rir(mixingTime:cutPoint);                % trim the tail
    %rirLate = rir(mixingTime:end); 
end

function [dB] = edc(signal)

    EDC = signal.^2; % get the output energy
    EDC = flip(EDC); % flip
    EDC = cumsum(EDC); % calculate cumulative sum 
    EDC = flip(EDC); % flip back
    EDC = EDC/max(abs(EDC)); % normalize
    EDC_dB = 10*log10((EDC)); % EDC in dB
    dB = EDC_dB;
end

