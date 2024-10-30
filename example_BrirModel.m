% Generate Binaural Dark Velvet Noise model of a target BRIR
% Jon Fagerström
% 5.6.2024

%% Init
clear; clc;
close all;
set(0,'defaultTextInterpreter','latex')
rng(420)
addpath('functions/');
addpath('classes/')
addpath('randpdf/')

%% PARAMS

fs = 48000;
tMix = 0.11;
tMixSamples = round(tMix*fs);

nFilter = 10;                          % number of dictionary filters
pulseFilterOrder = 5;                   % dictionary filter order   
preFilterOrder = 10;                    % pre filter order   

nWin = 2^12;                            % analysis window size
nHop = nWin/8;                          % analysis window hop size

density = 1500;                         % velvet-noise density (pulses/s)                
lenMod = 1;                             % broadband RT modifier

% coherence analysis params
winLen = 512;
win = hann(winLen);
noverlap = 0.75*winLen;
nfft = 4096;

%% LOAD BRIR
filePath = 'data/Promenadi_BRIR_Target.wav';
[brir, temp] = audioread(filePath);
brir = resample(brir,fs,temp); % use consistent sample rate


%% INIT DvnReverb
binauralReverb = DvnReverb(fs,2);

%% FIT DvnReverb
binauralReverb.fitModel(brir,tMix,nWin,nHop,nFilter,[preFilterOrder,pulseFilterOrder]);

%% SYNTHESIZE MODEL
binauralReverb.initModel(lenMod,density);
binauralReverb.prepare();

%% COHERENCE ANALYSIS
lateTarget = brir(tMixSamples:end,:);
lateModel = binauralReverb.ir(tMixSamples:end,:);
[cohTarget,~,f] = coher(lateTarget,win,noverlap,nfft,fs);
cohModel = coher(lateModel,win,noverlap,nfft,fs);


%% PLOTS
% SETTINGS
h = 300; w = 400;  

font = 13; lw = 2;

fig1 = figure('Renderer', 'painters', 'Position', [1000 310 w 0.5*h]);
semilogx(f,cohTarget,'LineWidth',lw); hold on;
semilogx(f,cohModel,'LineWidth',lw);
legend('Target', 'Model')