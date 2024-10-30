% Velve Noise Generator
% Author: Jon Fagerström
% Date: 28.9.2022
%   Arguments
%       <Fs>: Sampling rate
%       <density> Density [impulses/second]
%       <Ls> Sequence length [samples]
%   Returns
%       <k> pulse locations with shape (M,1)
%       <s> pulse signs with shape (M,1)


function [k, s, g, seq] = vn(Fs, density, Ls)
    Td = Fs/density;                        % average spacing between impulses = grid size
    M = round(Ls/Td);                       % Total number of impulses in the sequence
    r1 = rand(M,1);                         % random number sequence
    r2 = rand(M,1);                         % random number sequence
    
    s = 2 * round(r1) - 1;                  % signs
    k = ceil((0:M-1)'.*Td + r2.*(Td-1));    % pulse locations
    k(1) = 1;

    alpha = (-log(10.^(-60/20)))/Ls; % decay constant
    g = exp(-alpha*k);
    
    
    seq = zeros(Ls,1);
    seq(k) = s;
    if numel(seq)>Ls
        seq = seq(1:Ls); % make sure the sequence length is not exceeded
    end
end
