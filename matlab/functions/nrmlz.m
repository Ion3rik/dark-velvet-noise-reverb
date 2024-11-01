% Normalizing function
% Author: Jon Fagerström
% Date: 17.2.2021
% Arguments:
%           <input>: input signal to be normalized
%           <type>: normalization type: 'energy', 'peak'
%           <rem_ref>: reference rms level
% Returns:
%           <output>: normalized signal
function output = nrmlz(input, type, rms_ref)
    numChannel = size(input,2);
    L = size(input,1);
    output = zeros(L, numChannel);
    switch type
        case 'energy'
            for i = 1:numChannel
                rms = sqrt(sum(input(:,i).^2) / L);
                if rms == 0
                    output = zeros(size(input));
                    return
                end
                output(:,i) = input(:,i) * (rms_ref / rms);
            end
        case 'peak'
            for i = 1:numChannel
                output(:,i) = input(:,i) / max(abs(input(:,i))) * rms_ref;
            end
    end
end