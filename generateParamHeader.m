% Script to save DVN Reverb params to c++ header file to be loaded to the
% real-time plugin
% Jon Fagerström
% 28.8.2024


%% Init
clear; clc;
close all;
addpath('functions/');
addpath('classes/')
addpath('randpdf/')

%% PARAMS
rng(5)
fs = 48000;
tMix = 0.11;
tMixSamples = round(tMix*fs);

nFilter = 12;                          % number of dictionary filters
pulseFilterOrder = 2;                   % dictionary filter order   
preFilterOrder = 10;                    % pre filter order   

nWin = 2^12;                            % analysis window size
nHop = nWin/8;                          % analysis window hop size

density = 1000;                         % velvet-noise density (pulses/s)                
lenMod = 1;                             % broadband RT modifier

% coherence analysis params
winLen = 512;
win = hann(winLen);
noverlap = 0.75*winLen;
nfft = 4096;

roomNum = 'Room1';

%% LOAD BRIR
filePath = 'data/Promenadi_BRIR_Target.wav';
[brir, temp] = audioread(filePath);
brir = resample(brir,fs,temp); % use consistent sample rate
%brir = brir./max(abs(brir));

%% INIT DvnReverb
binauralReverb = DvnReverb(fs,2);

%% FIT DvnReverb
binauralReverb.fitModel(brir,tMix,nWin,nHop,nFilter,[preFilterOrder,pulseFilterOrder]);

%% SYNTHESIZE MODEL
binauralReverb.initModel(lenMod,density);
binauralReverb.prepare();

filterActivations = zeros(nFilter,1);
for m = 1:binauralReverb.M
    for q = 1:nFilter
        if(binauralReverb.pulseFilter(m) == q)
            filterActivations(q) = filterActivations(q)+1;
        end
    end
end


%% WRITE THE C++ HEADER FILE
path = '/Users/fagersj2/Documents/modern-real-time-audio/projects/DarkVelvetReverb/';

headerFilename = [path roomNum, '.h'];
fid = fopen(headerFilename, 'w'); % Open the file for writing

% Write the header guards and include directives
fprintf(fid, '#pragma once\n\n');
fprintf(fid, '#include <vector>\n');
fprintf(fid, '#include "JuceHeader.h"\n\n');



% Start writing the struct 
fprintf(fid, ['namespace Params\n{\n']);
fprintf(fid, ['struct ' roomNum 'Params\n{\n']);


% juce::AudioBuffer<unsigned int> pulseLocation;
% juce::AudioBuffer<float> pulseGain;
% std::vector<float> channelGain;
% std::vector<std::vector<float>> dictionaryFilterCoeff;
% std::vector<float> postFilterCoeff;

% Write the variables
M = binauralReverb.M;
nFilter = binauralReverb.nFilter;
nChannel = binauralReverb.nChannel;

fprintf(fid, 'unsigned int numFilter = %u; \n', nFilter);
fprintf(fid, 'unsigned int numPulse = %u; \n', M);
fprintf(fid, 'unsigned int numChannel = %u; \n', nChannel);


% filterRouting
fprintf(fid, 'std::vector<unsigned int> filterRouting = {');
for m = 1:M
    if m < M
        fprintf(fid, '%u, ', binauralReverb.pulseFilter(m)-1);  % Write each element with a comma
    else
        fprintf(fid, '%u', binauralReverb.pulseFilter(m)-1);    % Last element without comma
    end
end
fprintf(fid, '};\n\n');

% pulse locations
fprintf(fid, 'std::vector<std::vector<unsigned int>>  pulseLocation = {');
for ch = 1:nChannel
    fprintf(fid, '{');
    for m = 1:M
        if m < M
            fprintf(fid, '%u, ', binauralReverb.pulseLoc(m,ch)-1);  % Write each element with a comma
        else
            fprintf(fid, '%u', binauralReverb.pulseLoc(m,ch)-1);    % Last element without comma
        end
    end
    if ch < nChannel
        fprintf(fid, '},');
    else
        fprintf(fid, '}');
    end
end
fprintf(fid, '};\n\n');

% pulse gains
fprintf(fid, 'std::vector<std::vector<float>>  pulseGain = {');
for ch = 1:nChannel
    fprintf(fid, '{');
    for m = 1:M
        if m < M
            fprintf(fid, '%u, ', binauralReverb.pulseSign(m).*binauralReverb.pulseGain(m));  % Write each element with a comma
        else
            fprintf(fid, '%u', binauralReverb.pulseSign(m).*binauralReverb.pulseGain(m));    % Last element without comma
        end
    end
    if ch < nChannel
        fprintf(fid, '},');
    else
        fprintf(fid, '}');
    end
end
fprintf(fid, '};\n\n');

% channel gains
fprintf(fid, 'std::vector<float>  channelGain = {');
for ch = 1:nChannel
    if ch < nChannel
        fprintf(fid, '%u, ', binauralReverb.channelGain(ch));  % Write each element with a comma
    else
        fprintf(fid, '%u', binauralReverb.channelGain(ch));    % Last element without comma
    end
end
fprintf(fid, '};\n\n');

% dictionary filter coeffs
fprintf(fid, 'std::vector<std::vector<float>>  dictionaryFilterCoeff = {');
for q = 1:nFilter
    fprintf(fid, '{');
    for o = 1:pulseFilterOrder+1
        if o < pulseFilterOrder+1
            binauralReverb.filterCoeff.pulseA(o+1,q)
            fprintf(fid, '%u, ', binauralReverb.filterCoeff.pulseA(o+1,q));  % Write each element with a comma
        else
            binauralReverb.filterCoeff.pulseB(q)
            fprintf(fid, '%u', binauralReverb.filterCoeff.pulseB(q));    % Last element without comma
        end
    end
    if q < nFilter
        fprintf(fid, '},');
    else
        fprintf(fid, '}');
    end
end
fprintf(fid, '};\n\n');

% post filter coeffs
fprintf(fid, 'std::vector<float>  postFilterCoeff = {');
for o = 1:preFilterOrder+1
    if o < preFilterOrder+1
        fprintf(fid, '%u, ', binauralReverb.filterCoeff.preFilter(o+1));  % Write each element with a comma
    else
        fprintf(fid, '%u', 1);    % Last element without comma
    end
end
fprintf(fid, '};\n\n');

% Close the struct and namespace initialization 
fprintf(fid, '};\n\n');
fprintf(fid, '};\n\n');

% Close the file
fclose(fid);




%% PLOT
k = binauralReverb.pulseLoc;
gs = binauralReverb.pulseGain .* binauralReverb.pulseSign;
pulseSequence = zeros(max(max(k)),nChannel);

for ch = 1:nChannel
    pulseSequence(k(:,ch),ch) = gs(:,ch);
end

l = length(pulseSequence);
t = 1000*linspace(0,l/fs,l);

figure; 
plot(t,pulseSequence); xlim([0 371.52]); 


late = brir(tMixSamples:end,:);
late = late(1:l,:);

audiowrite('late.wav',late,fs);

%% POST FILTER TEST

f = linspace(0,fs/2,4096);

figure; 


%%
modi = 0;
for q = 1:nFilter
    B = binauralReverb.filterCoeff.pulseB(q);
    A = binauralReverb.filterCoeff.pulseA(:,q);
    A = A + modi;
    B = B;
    H = freqz(B, A,f, fs); 
    rms(impz(B,A))
    semilogx(f,db(abs(H))); hold on;
end



