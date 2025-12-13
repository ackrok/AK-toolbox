function [varargout] = phaseplot(LFP,st,Fs,binSize)
%phaseplot
%
%   Description: This function performs a Hilbert Transformation of a LFP
%   and uses that to find an LFP's instantaneous phase over time. Then it
%   utilizes a cell of multiple spike times or a single vector of spikes
%   times to create a phase distribution of firing rates. When no output
%   variables are specified the function creates a plot of the distribution
%
%   phaseplot(LFP,stCell,Fs,binSize) -- Plots the signal no output
%
%   [stPhaseDist_y,stPhasePlot_x] = phaseplot(LFP,stCell,Fs,binSize) -- If
%   you desire the functions output
%
%   Input:
%   LFP - A vector of the full length Local Field Potential - Should not be
%   downsampled
%   st - A cell array of multiple spikes times or a single vector of spike
%   times
%   Fs - Sampling Frequency of the acquisition - Should not be downsampled
%   sampling frequency
%   binSize - Bin Size for the phase distribution
%
%   Output:
%   stPhaseDist_y - The distribution of spikes within a phase
%   stPhaseDist_x - The x axis of the histogram
%
%   Author: Pratik Mistry, 2020

LFP_hilb = hilbert(LFP);
instPhase = (180/pi)*angle(LFP_hilb);
edgeVec = [-1*180-binSize:binSize:180+binSize];
stPhaseDist_y = zeros(1,length(edgeVec)-1);
stPhaseDist_x = edgeVec(1:end-1)+mean(diff(edgeVec))/2;
if ~iscell(st)
    tempSt = ceil(Fs*st);
    stPhase = instPhase(tempSt);
    stPhaseDist_y = histcounts(stPhase,edgeVec,'Normalization','probability');
else
    nSt = length(st);
    stPhaseDist_y = zeros(nSt,length(edgeVec)-1); 
    for n = 1:nSt
        tempSt = ceil(Fs*st{n});
        stPhase = instPhase(tempSt);
        %stPhaseDist_y = stPhaseDist_y + histcounts(stPhase,edgeVec);
        stPhaseDist_y(n,:) = histcounts(stPhase,edgeVec,'Normalization','probability');
    end
end
if nargout == 0
    figure; 
    if ~iscell(st)
        bar(stPhaseDist_x,stPhaseDist_y);
        xlim([-180 180]);
        xlabel('Phase'); ylabel('Probability');
    else
        plotM = floor(sqrt(nSt)); plotN = ceil(nSt/plotM); 
        for n = 1:nSt
            subplot(plotM,plotN,n);
            bar(stPhaseDist_x,stPhaseDist_y(n,:));
            xlim([-180 180]);
            xlabel('Phase'); ylabel('Probability');
        end
    end
else
    varargout{1} = stPhaseDist_y;
    varargout{2} = stPhaseDist_x;
end
end