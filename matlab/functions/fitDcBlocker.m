% Function for fitting a DC Blocker post filter onto a target response
% Jon Fagerström
% 7.6.2023

function [b,a] = fitDcBlocker(target, f, fs,order)
    % Analyse target response
    Target = abs(freqz(target, 1, f, fs))'; % magnitude response
    Target = octaveSmooth(Target,f,1); % octave smoothing
    Target = db(Target); % db
    Target = Target - max(Target); % normalize
    
    % DESIGN POST FILTER
    if order== 1
        refLevel = 3;
    elseif order == 2
        refLevel = 6;
    end
    [idx] = find((abs(max(Target)-Target)) < refLevel, 1, 'first'); % estimate xover frequency
    fc = f(idx);
    [b,a] = dcBlocker(fc,fs);
    
%     % DEBUG PLOT
%     PostFilter = db((freqz(b,a,f,fs)').^2);
%     figure; 
%     semilogx(f,Target, 'LineWidth',2); hold on;
%     semilogx(f, PostFilter, 'LineWidth',2);
%     xlim([20 20000]);
%     xlabel("Frequency (Hz)"); ylabel("Magnitude (dB)")
%     legend('Target', 'Post Filter')
end