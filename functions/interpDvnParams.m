% Function for interpolating analysis data to obtain DVN params

function [dataIntrp] = interpDvnParams(t, data, queryPoints)
    numFilters = size(data,1);
    M = numel(queryPoints);
    dataIntrp = zeros(numFilters,M);

    % Add extremity point to avoid extrapolation
    data = [data, data(:,end)];
    t = [t, queryPoints(end)];

    % Interpolate
    for n = 1:numFilters
        dataIntrp(n,:) = interp1(t,data(n,:),queryPoints, 'linear')'; % query point interpolation
    end
end