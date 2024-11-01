function [sLin, t, f] = plotSTFT(x, dynamicRange, frequencyRange, timeRange, fs, fftLength)
    windowLength = 2^9;
    win = hann(windowLength);
    overlap = windowLength /2;
    
    [sLin,f,t] = spectrogram(x,win,overlap,fftLength,fs); % compute spectrogram
    s = db(sLin); % convert to dB
    s = s - max(max(s)); % normalize to 0 dB
    s(s==-inf) = dynamicRange(1);
    sur = surf(t,f,s);
    sur.EdgeColor = 'none';
    axis xy; axis tight; view(0,90); % change to 2 D view

    %colormap(inferno()) % set color map
    cb = colorbar;
    set(cb,'YTick',-60:20:0,'YTicklabel',{'-60 dB', '-40 dB', '-20 dB', '0 dB'})
    if size(x,2) == 2
        set(cb,'YTick',0:5:20,'YTicklabel',{'0 dB', '5 dB', '10 dB', '15 dB', '20 dB'})
    end
    caxis(dynamicRange)
    set(gca,'TickLength',[0 0])

    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    xlim(timeRange)
    ylim(frequencyRange)
    ax = gca;
    ax.FontSize = 14;
    set(gca,'YScale', 'log', 'YTick',[50 250 1000 5000 20000], 'YTicklabel',{'50', '250' ,'1k', '5k', '20k'});
    grid off
end



