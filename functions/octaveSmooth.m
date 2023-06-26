
% FUNCTION FOR SMOOTHING FREQUENCY SPECTRUM
% Jon Fagerström
% 16.1.2023

function [spectrumSmoothed] = octaveSmooth(spectrum, f, N)

    [nBins, nSig] = size(spectrum);

    %% Smoothing Window
    fup = zeros(nBins,1); flo = zeros(nBins,1);
    for k = 1:length(f)
        fup(k,:) = f(k) * 2^(1/(2*N));
        flo(k,:) = f(k) / 2^(1/(2*N));
    end

    Win = zeros(nBins,nBins);
    for k = 1:nBins
        ix_lo   = find(f >= flo(k),1,'first');
        ix_up   = find(f <= fup(k),1,'last');
        Win(k,ix_lo:ix_up) = 1 / (ix_up - ix_lo + 1);
    end

%% Smooth
spectrumSmoothed = Win * spectrum ;

