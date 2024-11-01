function reflectionFilter = assignFiltersGreedy(p)
    % greedy assignment which should also work for changing probabilities
    numFilters = size(p,1); % number of filters
    numProbs = size(p,2);   % number of probabilities
    reflectionFilter = zeros(1,numProbs);

    % keep track of the last time a filter was chosen
    timeSinceLastPulse = 1000 * ones(numFilters,1);

    % go through the pulses sequentially
    for m = 1:numProbs
        % choose best filter for regularity my maximising the weighted
        % probability, i.e., the longer the time since last pulse the higher the chance. 
        r = rand(numFilters,1) * 0.1; % possibly some randomness
        [~,bestFilterIndex] = max((timeSinceLastPulse + r) .* p(:,m));
        reflectionFilter(m) = bestFilterIndex;
        
        % Reset timer and increment
        timeSinceLastPulse(bestFilterIndex) = 0;
        timeSinceLastPulse = timeSinceLastPulse + 1; 

        % Optionally, it's possible to cap the weighting
        % otherwise, if a previous filter which was at 0 probability turns
        % to non-zero, the filter would instantly be chosen. That can be
        % diserable or undesirable.
        % If capped too much, low probability filters will never play
        maxTime = 10000;
        timeSinceLastPulse(timeSinceLastPulse > maxTime) = maxTime;
    end


end