function  T60 = getT60(signal,tSTFT)
    nBins = size(signal,1);
    for k = 1:nBins
        
        EDC = abs(signal(k,:)).^2;          % get the output energy
        EDC = flip(EDC);                    % flip
        EDC = cumsum(EDC);                  % calculate cumulative sum 
        EDC = flip(EDC);                    % flip back
        EDC = EDC/max(abs(EDC));            % normalize
        EDC_dB = 10*log10((EDC));           % EDC in dB
        temp = find(EDC_dB<= -5, 1);        % in frames
        temp = find(EDC_dB(temp:end)<= -30, 1);
        T60(k) = 2*tSTFT(temp);             % in seconds
    end
    % Thanks Karolina Prawda!
    T60 = movmean(T60, 400);
end