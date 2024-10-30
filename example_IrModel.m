% Generate Extended Dark Velvet Noise model of a target IR
% Jon Fagerström
% 3.10.2023
%% Init
%clear; clc;
close all;
set(0,'defaultTextInterpreter','latex')
rng(420)
addpath('functions/');
%% Parameters
% GENERAL PARAMS
fs = 48000;                             % sample rate
% ANALYSIS PARAMS
winLen = 2^12;                          % analysis window size
hopSize = winLen/2;                     % analysis window hop size
fftLength = 4096;                       % fft length for the analysis
fRange = [20, 20000];                   % frequency range for the analysis                         
% SYNTHEIS PARAMS
nFilters = 10;                          % number of dictionary filters
pulseFilterOrder = 2;                   % dictionary filter order   
preFilterOrder = 10;                    % pre filter order     
density = [1500, 1500];                  % velvet-noise density as a range [startDensity, endDensity], (pulses/s)                
lenMod = 1;                             % broadband RT modifier
% PLOT PARAMS
dynamicRange = [-60, 0];                % spectrogram dynamic range
frequencyRange = [50 20000];            % spectrogram frequency range

%% LOAD TARGET RESPONSE
filePath = 'soundExamples/Promenadi Hall/Promenadi_IR_Target.wav';
[rirTarget, temp] = audioread(filePath);
rirTarget = resample(rirTarget,fs,temp); % use consistent sample rate

%% ESTIMATE MIXING TIME
tMix = 1;  % mixing time in samples

%% EXTRACT MONO LATE REVERB

[lateTarget, earlyTarget] = extractLateReverb(rirTarget(:,1),...
                                              tMix,...
                                              80); 
lateTargetOg = lateTarget;                      % save original late target
timeRange = [0 1] .* length(lateTarget)/fs;     % update plot params    
lenLate = length(lateTarget);                   % length of late reverberation                               
lenEarly = numel(earlyTarget);                  % length of direct+early 
lenRIR = lenEarly + lenLate;                    % total length of target IR

%% PRE-FILTER
preFilter = lpc(lateTarget(1:winLen) .* hann(winLen),preFilterOrder);   % estimate pre filter from the first frame
lateTarget = filter(preFilter,1,lateTarget);                            % apply pre filter

%% ANALYSIS 
win = hann(winLen); % analysis window
[targetResponses,f,tFrame, lp] = analyzeRir(lateTarget,... 
                                            fs,...
                                            win,...
                                            hopSize,...
                                            fftLength,...
                                            pulseFilterOrder);  % run analysis
nFrames = size(targetResponses,2);                              % number of STFT frames
nBins = numel(f);                                               % number of frequency bins

%% DESIGN PULSE FILTERS

logSubSet = unique(round(logspace(log10(1),...
                   log10(nFrames),...
                   nFilters)));                 % pick subset of analysed filters, spaced logartihmically in time
A = lp(logSubSet,:)'; 
B = ones(1,nFilters);
nFilters = size(A,2);                           % number of dictionary filters

filterResponses = zeros(nBins,nFilters);
for n = 1:nFilters
    ir = impz(B(:,n), A(:,n));
    enrg = sqrt(sum(ir.^2));
    B(:,n) = B(:,n) ./ enrg;                    % normalize filter energy
    filterResponses(:,n) = abs(freqz(B(:,n),...
                                     A(:,n),...
                                     f,...
                                     fs))';     % compute magnitude spectra
end

%% ESTIMATE FILTER PROBABILITIES

[x, xSparse, fittedResponses, gNormFrame] = estimateFilterProbabilities(targetResponses,...
                                                                        filterResponses);       % compute filter activations
p = x./gNormFrame;                                                                              % convert to probability

%% SYNTHESIS
lenModel = round(lenLate*lenMod);
tFrame = round(tFrame * lenMod);

% GENERATE VELVET NOISE %
[k,s] = vnTimeVar(fs,density,lenModel);                         % generate pulse locations and signs
M = length(k);                                                  % number of pulses in the sequence

% ASSIGN PULSE FITERS %
pIntrp = interpDvnParams(tFrame,p,k);                           % interpolate probabilities
filtList = assignFiltersGreedy(pIntrp);                         % assign filters based on the probabilities

% ESTIMATE PULSE GAINS %
gNormPulse = interpDvnParams(tFrame,...
                             gNormFrame,...
                             k);                                % interpolate gains
g = gNormPulse .* s;                                            % combine with pulse signs

% COMPUTE IMPULSE RESPONSE %
[lateModel, dvnFilterOut] = convDvn(1,...
                                    k,...
                                    g,...
                                    filtList,...
                                    B,...
                                    A);                         % compute model impulse response 

lateModel = [lateModel; zeros(lenModel-length(lateModel),1)];   % force same length as target

%% POST FILTER   
lateModel = filter(1,preFilter,lateModel);      % apply inverse of the pre filter
[b,a] = fitDcBlocker(lateTargetOg, f, fs, 1);   % fit the dc blocker
lateModel = filter(b,a,lateModel);              % apply dc bolcker


%% CONSTURCT FULL IR
eNorm = rms(lateTargetOg);                      % late target rms level
lateModel = nrmlz(lateModel, 'energy', eNorm);  % normalize energy to match target
model = [earlyTarget; lateModel];               % construc full mode IR
target = [earlyTarget; lateTargetOg];           % reconstruct full target IR

%% PLOTTING

% SETTINGS
h = 300; w = 400;
dbSpread = 100;     
colorMod = linspace(1.2,0.4,nFilters);
c = [0, 0.4470, 0.7410];
cEarly = [.5 .5 .5];
cFDN = [0.9290 0.6940 0.1250];
colorMap = [c(1) * colorMod; c(2)*colorMod; c(3)*colorMod];
font = 13; lw = 2;

% PRE PROCESS
model_dB = db(model);
target_dB = db(target);
earlyTarget_dB = db(earlyTarget); refLevel = max(target_dB); earlyTarget_dB = earlyTarget_dB - refLevel;
lateModel_dB = db(lateModel); lateModel_dB = lateModel_dB - refLevel;
lateTarget_dB = db(lateTargetOg); lateTarget_dB = lateTarget_dB - refLevel;


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
tModelLate = linspace(tMix/fs,lenLate/fs,lenModel);

plot(tEarly,earlyTarget_dB, 'color', cEarly); hold on;
plot(tModelLate, lateModel_dB, 'k'); tSettings(font,'Time (s)', 'Magnitude (dB)');
xlim([-0.05 lenLate/fs]); 
ylim([-60 0]);




%% FUNCTIONS

function fSettings(font)
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    xlim([20 20000])
    ax = gca;
    ax.FontSize = font;
    set(gca,'XScale', 'log', 'XTick',[ 20 50 100 250 500 1000 2000 4000 8000 16000], 'XTicklabel',{'20','50', '100', '250', '500', '1k', '2k', '4k', '8k', '16k'});
    %grid on;

end

function plotProbMap(p)
    imagesc((p));
    xlabel('Time (frames)')
    ylabel('Filter')
    ax = gca;
    ax.FontSize = 12;
    cb = colorbar;
    cb.Label.String = 'Probability';
    cb.Label.Interpreter = 'latex';
    set(cb,'FontSize',12);
    load("colorMapOranges.mat",'oranges');
    colormap(oranges) % set color map
    set(gca, 'clim', [0 1]); 
    set(gca,'YDir','normal')
end

function tSettings(font, x_label, y_label)
    xlabel(x_label)
    ylabel(y_label)
    ax = gca;
    ax.FontSize = font;
end

function  T60 = getT60(signal,tSTFT)
    nBins = size(signal,1);
    for k = 1:nBins
        
        EDC = abs(signal(k,:)).^2;          % get the output energy
        EDC = flip(EDC);                    % flip
        EDC = cumsum(EDC);                  % calculate cumulative sum 
        EDC = flip(EDC);                    % flip back
        EDC = EDC/max(abs(EDC));            % normalize
        EDC_dB = 10*log10((EDC));           % EDC in dB
        temp = find(EDC_dB<= -5, 1);        % in frames
        temp = find(EDC_dB(temp:end)<= -30, 1);
        T60(k) = 2*tSTFT(temp);             % in seconds
    end
    % Thanks Karolina Prawda!
    T60 = movmean(T60, 400);
end

