% Extract time dependant magnitude responses of the RIR
% Jon Fagerström
% Updated 28.9.2022

function [s,f,t, lp, gDenoise] = analyzeRir(rir, fs, win, hopSize, NFFT, lpOrder)
    %% PARAMS & INIT
    L = length(rir);                        % length of the analysed RIR in samples
    nWin = length(win);                     % window length
    %lpOrder = 10;                           % LP order for smoothing
    f = linspace(0,fs/2,NFFT/2+1);           % frequency axis
    nFrames = floor((L-nWin)/hopSize+1);    % number of frames
    s = zeros(NFFT/2+1,nFrames);               % init spectrogram matrix
    t = zeros(1,nFrames);                   % init time vector (samples)

    pin = 0;                                % counter for sliding the frame
    rir = [rir;zeros(nWin/2,1)];            % zero pad the input
    lp = zeros(nFrames,lpOrder+1);          % init lp filter matrix
    %% DENOISING
%     if denoise
%         noiseEstimate = rir(end-4800+1:end); % estimate noise from the tail
%         NoiseEstimate = (abs(freqz(noiseEstimate,1,f,fs)))'; 
%         %NoiseEstimate = octaveSmooth(NoiseEstimate, f', 3);
%     end
    %% Compute STFT
    for n = 1:nFrames
        grain = rir(pin+1:pin+nWin).*win;         % extract frame and apply window function
        lp(n,:) = lpc(grain, lpOrder);            % compute LP
        %s(:,n) = (abs(freqz(1,lp(n,:),f,fs)))';  % compute mag response
        temp = fft(grain, NFFT);                        % compute FFT
        s(:,n) = abs(temp(1:NFFT/2+1));           % compute one-sided mag response           
        %s(:,n) = octaveSmooth(s(:,n), f', 3);           % octave smoothing
        t(n) = pin;               % save frame time
        pin = pin + hopSize;                      % make the hop
        %norm(grain)
    end
    %t(1) = 1; t(end) = L;
    gDenoise = [];
%     if denoise
%         figure; 
%         subplot(3,1,1); semilogx(f,db(NoiseEstimate)); title('Noise estimate')
%         subplot(3,1,2); semilogx(f,db(s)); title('Rir Spectra')
%         sOg = s;
%         s = max(0,s - NoiseEstimate); % remove noise spectrum
%         gDenoise = rms(sOg.^2) - rms(s.^2); % energy compensation for the denoising
%     
%         subplot(3,1,3); semilogx(f,db(s)); title('Rir Spectra (Denoised)')
%     end
end

