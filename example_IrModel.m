% Generate Extended Dark Velvet Noise model of a target late-reveberation IR
% Jon Fagerström
% 3.10.2023

clear; clc;
close all;

%% PARAMETERS

% GENERAL PARAMS
fs = 48000;                             % sample rate

% ANALYSIS PARAMS
nWin = 2^12;                            % analysis window size
nHop = nWin/2;                          % analysis window hop size
fftLength = 4096;                       % fft length for the analysis
fRange = [20, 20000];                   % frequency range for the analysis 

% SYNTHESIS PARAMS
nFilter = 10;                           % number of dictionary filters
pulseFilterOrder = 2;                   % dictionary filter order   
preFilterOrder = 10;                    % pre filter order     
density = 1500;                         % velvet-noise density (pulses/s) 

% PLOT PARAMS
dynamicRange = [-60, 0];                % spectrogram dynamic range
frequencyRange = [50 20000];            % spectrogram frequency range

%% LOAD TARGET RESPONSE
filePath = 'soundExamples/Promenadi Hall/Promenadi_IR_Target.wav';
[target, temp] = audioread(filePath);
target = resample(target,fs,temp); % use consistent sample rate

%% ESTIMATE MIXING TIME
tMix = round(0.11*fs);  % mixing time in samples

%% INIT DvnReverb OBJECT
monoReverb = DvnReverb(fs,1);

%% FIT TO TARGET RIR
monoReverb.fitModel(target,tMix,nWin,nHop,nFilter,[preFilterOrder,pulseFilterOrder]);

%% SYNTHESIZE MODEL
monoReverb.initModel(1,density);
monoReverb.prepare()

%% SAVE IMPULSE RESPONSES
% model
lateModel = monoReverb.lateModel;
model = monoReverb.ir; 

% target
lateTarget = monoReverb.lateTarget;
earlyTarget = monoReverb.earlyTarget;
target = [earlyTarget; lateTarget];

%% PLOTTING

% SETTINGS
h = 300; w = 400;
dbSpread = 100;     
colorMod = linspace(1.2,0.4,nFilter);
c = [0, 0.4470, 0.7410];
cEarly = [.5 .5 .5];
cFDN = [0.9290 0.6940 0.1250];
colorMap = [c(1) * colorMod; c(2)*colorMod; c(3)*colorMod];
font = 13; lw = 2;

timeRange = [0, length(model)/fs];
lenLate = length(lateModel);
lenEarly = length(earlyTarget);

% PRE PROCESS
model_dB = db(model);
target_dB = db(target);
earlyTarget_dB = db(earlyTarget); refLevel = max(target_dB); earlyTarget_dB = earlyTarget_dB - refLevel;
lateModel_dB = db(lateModel); lateModel_dB = lateModel_dB - refLevel;
lateTarget_dB = db(lateTarget); lateTarget_dB = lateTarget_dB - refLevel;


% Plot Target spectrogram
fig1 = figure('Position', [400 1 w 0.525*h]);
[sTarget, tSTFT, fSTFT] = plotSTFT(target, dynamicRange, frequencyRange, timeRange, fs, fftLength); hold on;
T60Target = getT60(sTarget, tSTFT); semilogy(T60Target,fSTFT,'--w', 'linewidth', 2)

% Plot Model spectrogram
fig2 = figure('Position', [400 300 w 0.525*h]);
sModel = plotSTFT(model, dynamicRange, frequencyRange, timeRange, fs, fftLength); hold on;
T60Model = getT60(sModel, tSTFT);
semilogy(T60Model,fSTFT, 'w', 'linewidth', 2);
semilogy(T60Target,fSTFT, '--w', 'linewidth', 2);

% Plot Target IR
irLim = max(abs([target;model;]));
tTargetLate = linspace(tMix/fs,lenLate/fs,lenLate);
tEarly = linspace(0,lenEarly/fs,lenEarly);
fig3 = figure('Renderer', 'painters', 'Position', [1000 310 w 0.5*h]);
plot(tEarly, earlyTarget_dB, 'color', cEarly); hold on;
plot(tTargetLate, lateTarget_dB, 'k'); tSettings(font,'Time (s)', 'Magnitude (dB)');
xlim([-0.05 lenLate/fs]); %ylim([-irLim, irLim]);
ylim([-60 0]);

% Plot Model IR
fig4 = figure('Renderer', 'painters', 'Position', [1000 610 w 0.5*h]);
tModelLate = linspace(tMix/fs,lenLate/fs,lenLate);

plot(tEarly,earlyTarget_dB, 'color', cEarly); hold on;
plot(tModelLate, lateModel_dB, 'k'); tSettings(font,'Time (s)', 'Magnitude (dB)');
xlim([-0.05 lenLate/fs]); 
ylim([-60 0]);




%% FUNCTIONS

function tSettings(font, x_label, y_label)
    xlabel(x_label)
    ylabel(y_label)
    ax = gca;
    ax.FontSize = font;
end



