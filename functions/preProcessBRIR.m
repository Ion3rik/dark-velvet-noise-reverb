% Pre-process BRIR
% Remove silent beginning and cut to dynamic range

function [early,late] = preProcessBRIR(brir,mixingTime,dynamicRange,trim)
    late = brir;
    nChannel = size(brir,2);
    if trim
        startIdx = zeros(nChannel,1);
        for ch = 1:nChannel
            idx = find(late(:,ch).^2 > 0.001,1); % find where the signal energy is more -60 dB
            if isempty(idx)
                startIdx(ch) = 1;
            else
                startIdx(ch) = idx;         
            end
        end
            startIdx = min(startIdx); % pick the channel that starts earlier
            late = late(startIdx:end,:);                        % remove silent begining
    end
    early = late(1:mixingTime-1,:);
    if ~isempty(dynamicRange)
        brirEdc = edc(late);                                % compute the energy decay
        cutPoint = find(brirEdc <= -dynamicRange,1);        % find where to cut the tail based on the dynamic range
        late = late(mixingTime:cutPoint,:);                   % trim the tail
    else
        late = late(mixingTime:end,:);                   % trim the tail
    end
end

function [dB] = edc(signal)

    EDC = sum(signal.^2,2); % get the output energy
    EDC = flip(EDC); % flip
    EDC = cumsum(EDC); % calculate cumulative sum 
    EDC = flip(EDC); % flip back
    EDC = EDC/max(abs(EDC)); % normalize
    EDC_dB = 10*log10((EDC)); % EDC in dB
    dB = EDC_dB;
end