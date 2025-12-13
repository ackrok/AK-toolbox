function fig = plotUnitProp(data) 
%Plot waveform and spiking properties, to be used for unit characterization
%
% fig = plotUnitProp(data)
%
% Description: This function plots waveform half-width and peak-latency,
% unit firing rate, and proportion of recording when ISI > 2s for all units
%
% INPUT
%   'data' - structure with sub-structures data.WF and data.spk
%       must have previously run processWF, getWFprop, getSpikeProp
%
% OUTPUT
%   'fig' - optional return figure hanle
%
% Created by: Anya Krok, March 2020

    if ~isfield(data.WF,'halfWidth')
        error('Waveforms are not processed. Run processWF and getWFprop.');
    end
    
    halfWidth   = [data.WF.halfWidth]; 
    peakLatency = [data.WF.peakLatency]; 
    amplitude   = [data.WF.amplitude]; 
    fr          = [data.spk.fr]; 
    pISI2       = [data.spk.pISI2];
    
    nUnits = length(data.clusters);
    plotme = find(~isnan(halfWidth)); 

    fig = figure; fig.Position(3) = 1650;
    c = jet(nUnits);
    
    subplot(2,3,[1,4]);
        scatter(halfWidth,peakLatency,[],c,'filled');
        title(sprintf('%s-%s - Waveform Properties',data.mouse,data.date))
        xlabel('half-width (ms)'); ylabel('peak latency (ms)'); 
        xlim([0 4]); ylim([0 4])

    subplot(2,3,2); 
        b1 = bar(halfWidth,1); b1.FaceColor = 'flat';
        title('Waveform Half Width'); 
        ylabel('half-width (ms)'); %xlabel('unit number'); 

    subplot(2,3,5); 
        b2 = bar(peakLatency,1); b2.FaceColor = 'flat'; 
        title('Waveform Peak Latency');  
        ylabel('peak latency (ms)'); xlabel('unit number');

    subplot(2,3,[3,6])
        scatter(fr(plotme), pISI2(plotme),[],c(plotme,:),'filled');
        title('Firing Rate vs. Prop ISI > 2s')
        xlabel('firing rate (Hz)'); 
        ylabel('proportion of recording where ISI > 2s'); ylim([0 1])

    for x = 1:nUnits
        b1.CData(x,:) = c(x,:); b2.CData(x,:) = c(x,:);
    end
    
end