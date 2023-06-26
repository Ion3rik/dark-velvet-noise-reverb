% First order DC Blocker IIR
% Jon Fagerström
% 7.6.2023
% REF: FILTER-BASED ALIAS REDUCTION FOR DIGITAL CLASSICAL WAVEFORM
% SYNTHESIS bu Jussi Pekonen and Vesa Välimäki 
% https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4517564

function [b,a] = dcBlocker(fc, fs)
%     fc = 50;
%     fs = 48000;
    p = tan(pi/4 - pi*fc/fs);
    b = [1+p, -(1+p)];
    a = [2, -2*p];
    %f = linspace(20,20000,4096);
    %figure; semilogx(f,db(freqz(b,a,f,fs)'));
end