% DVN REVERBERATOR
% Author: Jon Fagerström 
% Date: 9.4.2024 
% Desicription: This class implements the Dark Velvet Noise Reverberator
% presented in JAES and DAFX24 (binaural version)

classdef DvnReverb < handle
    properties
        % GENERAL PARAMS %
        fs = 48000;                                 % sampling frequency
        density = 1500;                             % pulse density
        nFilter = 10;                              % number of filters in the DVN structure'
        nChannel = 1;                             % number of channels in the target rir
        mixingTime = 0;                             % point where to change between early reflection and late reverb processing
        nfft = 4096;                                % FFT length
        freqFft = linspace(0, 48000/2, 4096);   % vector of frequencies at which to run the analysis and fitting
        
        early;                                  % early part of the target RIR

        % FIT PARAMS
        prob;           % DVN filter probabilities per frame
        frameGain;      % Decay gains per frame
        frameTime;      % Analysis frame times
        nTarget;        % length of the target response
        targetRIR;      % target RIR

        % SYNTHESIS PARAMS %
        M;           % DVN number of pulses
        pulseGain;   % DVN pulse gains [M, 1]
        pulseFilter  % DVN pulse filter routing
        pulseLoc;    % DVN pulse locations
        pulseSign    % DVN pulse signs
        filterCoeff = struct(); % DVN filter coefficients (B,A,preFilter) 
        channelEnergy; % Left and right channel energies
        channelGain;

        % PROCESS MEMORY
        ir

    end
    
    methods
        %% CONSTRUCTOR
        function obj = DvnReverb(fs,nChannel)
            if nChannel > 2
                error("Supporting up to 2 channels, for now.")
            end
            obj.fs = fs;
            obj.freqFft = linspace(0, fs/2, obj.nfft); 
            obj.nChannel = nChannel;
        end
        %% FIT MODEL
        function fitModel(obj,rir, tMix, nWin, nHop,nFilter, filterOrder)
            if obj.nChannel < size(rir,2) && size(rir,2) < 3
                error("The target RIR has more channels than the model. Increase the model channels.")
            elseif obj.nChannel > size(rir,2)
                warning("Preparing a binauralized version of the mono RIR.")
            elseif size(rir,2) > 2
                error("Target RIR has too many channels. Supporting up to 2 channels, for now.")
            end
            obj.nFilter = nFilter;
            obj.targetRIR = rir; % save the target RIR

            [obj.early, rir] = preProcessBRIR(rir,round(obj.fs*tMix),60,true); % extract early and late parts
            obj.channelEnergy = rms(rir); % save energy 
            rir = rir(:,1); % use the left channel for filter extraction
            obj.nTarget = numel(rir); % length of the target rir
            win = hann(nWin); % analysis window
            
            % PRE-FILTER
            preFilter = lpc(rir(1:nWin) .* win, filterOrder(1));   
            rir = filter(preFilter,1,rir); 
            
            % ANALYSIS 
            [targetResponses,f,obj.frameTime, lp] = analyzeRir(rir, obj.fs, win, nHop, obj.nfft,filterOrder(2));
            nFrame = size(targetResponses,2);                              % number of STFT frames
            [b,a] = fitDcBlocker(rir, f, obj.fs, 2);      % fit the dc blocker

            % DESIGN PULSE FILTERS
            logSubSet = unique(round(logspace(log10(1), log10(nFrame), obj.nFilter)));
            A = lp(logSubSet,:)'; 
            B = ones(1,obj.nFilter);
            obj.nFilter = size(A,2); % number of filters might change slightly                          
            
            % COMPUTE FILTER RESPONSES
            filterResponses = zeros(numel(f),obj.nFilter);
            for n = 1:obj.nFilter
                ir = impz(B(:,n), A(:,n));
                enrg = sqrt(sum(ir.^2));
                B(:,n) = B(:,n) ./ enrg; % normalize energy                       
                filterResponses(:,n) = abs(freqz(B(:,n), A(:,n), f, obj.fs))'; % compute magnitude response    
            end
            % SAVE FILTER COEFFICIENTS
            obj.filterCoeff.pulseA = A;
            obj.filterCoeff.pulseB = B;
            obj.filterCoeff.dcBlockA = a;
            obj.filterCoeff.dcBlockB = b;
            obj.filterCoeff.preFilter = preFilter;

            % ESTIMATE FILTER PROBABILITIES
            [x, ~, ~, obj.frameGain] = estimateFilterProbabilities(targetResponses,filterResponses); 
            obj.prob = x./obj.frameGain;                                                                             
           
        end
       %% SET MODEL LENGTH
       function setLength(obj,lenMod)
            % APPLY THE LENGTH MODIFIER
            obj.nTarget = round(obj.nTarget * lenMod);
            obj.frameTime = round(obj.frameTime * lenMod);
       end
       %% SET PULSES
       function setPulses(obj, density)
            % COMPUTE PULSE LOCATIONS AND SIGNS
            Td = obj.fs/density;               
            obj.M = round(obj.nTarget/Td);
            maxJitter = round(0.001*obj.fs); % 1 ms is good starting point, change if needed

            if size(obj.targetRIR,2) == 2 && obj.nChannel == 2 % Target RIR has 2 channels, model has 2 channels
                jitter = estimateJitter(obj.targetRIR, maxJitter, obj.fs,[obj.M,1],0);
                [~,obj.pulseLoc,obj.pulseSign] = velvetJitter(obj.nTarget, Td, 0, jitter, "Additive");

            elseif size(obj.targetRIR,2) == 1 && obj.nChannel == 2 % Target RIR has 1 channel, model has 2 channels
                jitter = estimateJitter("theoretical", maxJitter, obj.fs,[obj.M,1],0); % use binaural diffuse field coherence
                [~,obj.pulseLoc,obj.pulseSign] = velvetJitter(obj.nTarget, Td, 0, jitter, "Additive");

            elseif size(obj.targetRIR,2) == 2 && obj.nChannel == 1 % Target RIR has 2 channel, model has 1 channels
                jitter = estimateJitter("theoretical", maxJitter, obj.fs,[obj.M,1],0);
                [~,obj.pulseLoc,obj.pulseSign] = velvetJitter(obj.nTarget, Td, 0, jitter, "Additive");
            elseif size(obj.targetRIR,2) == 1 && obj.nChannel == 1 % Target RIR has 1 channel, model has 1 channels
                [obj.pulseLoc,obj.pulseSign] = vn(obj.fs,density,obj.nTarget);
            end
            %obj.pulseSign = obj.pulseSign * sqrt(Td); % normalize with the density
       end
       %% SET PULSE GAINS
       function setPulseGains(obj)
            % ESTIMATE PULSE GAINS %
            obj.pulseGain = interpDvnParams(obj.frameTime, obj.frameGain,sort(round(mean(obj.pulseLoc,2)))); % interpolate         
       end

       %% SET PULSE FILTERS
       function setPulseFilters(obj)
            % ASSIGN PULSE FITERS %
            probIntrp = interpDvnParams(obj.frameTime,obj.prob,sort(round(mean(obj.pulseLoc,2)))); % interpolate
            obj.pulseFilter = assignFiltersGreedy(probIntrp); % obtain pulse filter routing                           
       end

       %% INIT MODEL
       function initModel(obj, lenMod, density) % calls the 4 set functions
           obj.setLength(lenMod);
           obj.setPulses(density);
           obj.setPulseGains();
           obj.setPulseFilters();
       end

       %% PREPARE
       function prepare(obj)
            
            obj.ir = zeros(obj.nTarget,obj.nChannel);
            obj.pulseGain = obj.pulseGain';
            for q = 1:obj.nFilter % for each dictionary filter
                idx = obj.pulseFilter == q;
                vq = zeros(obj.nTarget,obj.nChannel); % white sub sequence
                for ch = 1:obj.nChannel
                    vq(obj.pulseLoc(idx,ch),ch) = obj.pulseSign(idx,ch) .* obj.pulseGain(idx); % get white sub sequence
                    obj.ir(:,ch) = obj.ir(:,ch) + filter(obj.filterCoeff.pulseB(q), obj.filterCoeff.pulseA(:,q),vq(:,ch)); % add to the ir
                end
            end
            obj.ir = filter(1, obj.filterCoeff.preFilter,obj.ir); % post filter
            obj.ir = filter(obj.filterCoeff.dcBlockB, obj.filterCoeff.dcBlockA, obj.ir); % dc blocker
            obj.ir = obj.channelEnergy(ch) ./ rms(obj.ir) .* obj.ir; % apply channel gain
            obj.channelGain = obj.channelEnergy(ch) ./ rms(obj.ir);


            % construct full RIR
            obj.ir = [obj.early; obj.ir]; 
       end

       %% PROCESS
       function y = process(obj,x)
            y = conv2(x,obj.ir);
       end
    end

end