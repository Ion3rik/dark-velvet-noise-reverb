function [coh,cross,f] = coher(brir,win,noverlap,nfft,fs)
    HL = stft(brir(:,1),fs,"Window",win,"OverlapLength",noverlap,"FFTLength",nfft,"FrequencyRange","onesided");
    HR = stft(brir(:,2),fs,"Window",win,"OverlapLength",noverlap,"FFTLength",nfft,"FrequencyRange","onesided");
    
    SL = sum(abs(HL).^2,2); % PSD Left
    SR = sum(abs(HR).^2,2); % PSD Right

    cross = sum(HL.*conj(HR),2); % CPSD

    f = linspace(0,fs/2,length(cross));

    coh = abs(octaveSmooth(cross,f,3)).^2 ./ octaveSmooth(SL .* SR,f,3);
    %coh = abs(cross).^2 ./ (SL .* SR); % unsmoothed

end
