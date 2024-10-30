% Function for estimating corresponding jitter from a target BRIR
% Jon Fagerström
% 20.2.2024
% Arguments:
%           <targetBRIR> BRIR with the target coherence
%           <maxJitter> Maximum allowed jitter in samples
%           <dim>       Dimensions of the output
% Outputs:
%           <jitter>    matrix of jitter samples with dimensions specified
%           in <dim>

function jitter = estimateJitter(targetBRIR, maxJitter,fs,dim,plots)
    
    if strcmp(targetBRIR, 'theoretical')
        % replace target brir with the empirical diffuse noise
        load("data/binauralDiffuseNoise.mat","wnBinaural");
        targetBRIR = wnBinaural;

    elseif strcmp(targetBRIR, 'empirical')
        load('empirical_itds','itd');
        itd = itd(randi(length(itd),dim))'*1.0000e-06;
        jitter = round(itd*fs);
        return;
    end

    % % extend with noise
    %noise = rand(60*48000,1);
    %targetBRIR = conv2(noise,targetBRIR);
    % estimate cross power spectral density
    
    winLen = 512;
    noverlap = winLen*0.75;
    win = hann(winLen);
    nfft = 4096;

    for k = 1:size(targetBRIR,3)
        [~,crossCurrent] = coher(targetBRIR, win, noverlap, nfft, fs); % use coherence instead
        %cross = cpsd(targetBRIR(:,1), targetBRIR(:,2), win, noverlap, nfft); % use cpsd
        %nfft = 2*48000;
        %cross = fft(targetBRIR(:,1),nfft) .*conj(fft(targetBRIR(:,2),nfft));       % CPSD between left and right channel
    
        % octave smoothing
        %f = linspace(0,48000/2,numel(cross));
        %cross = smoothSpec(cross,f, 3); % apply third-octave smoothing
    
        % mooving average
        %cross = complex(movmean(real(cross),1000), imag(cross));
    
       
    cross(:,k) = crossCurrent;
    end
    cross = mean(cross,2);
    % make two sided
    cross = [cross; flip(conj(cross(2:end)))]; 

    % compute jitter distribution
    currentDist = fftshift(ifft(cross));                                  
    lenJitterDist = length(currentDist);                                                   
    currentDist = max(0,currentDist(round(lenJitterDist/2)-maxJitter:round(lenJitterDist/2)+maxJitter));   
    t = linspace(-maxJitter,maxJitter,length(currentDist));
    
    % smooth distribution
    %jitterDist = movmean(jitterDist,50);
    jitterDist = mean(currentDist,2);
    % draw samples from the distribution  
    jitter = (round(randpdf((jitterDist),t,dim))); 
    
    
    if plots
        figure;
        subplot(2,1,1)
        plot(t, currentDist)
        title('Designed Distribution'); 
        subplot(2,1,2)
        histNorm(jitter, 1000);
        title('Drawn Distribution')
    end
end

function itd = computeItd(phi, theta, freeField)

a = 0.09; % head radius in m
c = 343; % speed of sound in m/s

    if freeField % free field delays
        receivers(1,:) = [a, 0, 0];
        receivers(2,:) = [-a, 0, 0];
    
        dist(1,:) = sqrt(sum((phi - receivers(1,:)).^2,2));
        dist(2,:) = sqrt(sum((phi - receivers(2,:)).^2,2));
    
        itd = abs(dist(2,:) - dist(1,:)) ./ c;
        itd = itd';
        return;
    end

% woodworth approximation formula for itd
%itd = a / c * (pi - phi + sin(phi)) .* cos(theta);

%itd(phi < pi / 2)= (a / c * (phi(phi < pi / 2) + sin(phi(phi < pi / 2)))) .* cos(theta(phi < pi / 2));

phi = phi + pi / 2;

itd = a / c * (cos( phi ) + pi / 2 - phi ) .* cos(theta);

end