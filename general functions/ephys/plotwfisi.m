function [fig] = plotwfisi(data, varargin)
% plotUnits - Plot unit waveform and log-ISI distribution
%
%   [fig] = plotwfisi (data);
%   [fig] = plotwfisi (data, units);
%
%   Description: This function plots histogram of spike times over the
%   length of silicon probe recording along with average waveform and
%   log-ISI distribution for specified units. When available, velocity
%   signal and photometry signal will also be plotted.
%
%   Input:
%   - data - Structure containing informatin about silicon probe recording
%         data.WF: unit waveforms
%         data.spk: unit temporal properties, including ISI
%   - units (optional) - Vector containing unit numbers to plot
%   
%   Output:
%   - fig - Figure handle 
%
% Author: Anya Krok, 2019
%

%% 
    if ~isfield(data,'WF')
        error('Error: no waveforms in data structure.')
        return
    end
    if ~isfield(data,'spk')
        error('Error: no ISI in data structure.')
        return
    end
    
    if nargin == 2
        plotme = varargin{1}; %plot only specified units
    else
        plotme = 1:length(data.clusters); %plot all units in recording
    end

    fig = figure; spN = 2; spM = length(plotme); %set figure properties
    
%% Plot
    for x = 1 : length(plotme)
        n = plotme(x);
        
        spwf(x) = subplot(spM, spN, (x*spN)-1);
            [z_wf, mu] = zscore(data.WF(n).topMu); 
            z_wfsem = data.WF(n).topSem./mu; %pull from data structure
            shadederrbar(data.WF(n).time, z_wf, z_wfsem, 'k'); %plot average and stdev of waveform from max channel
            title(sprintf('unit #%d - WF',n));  ylabel('Amp (SD)');
            xlim([-3 3])

        spisi(x) = subplot(spM, spN, (x*spN));
            unitISI = data.spk(n).ISI; %pull from data structure
            [~, edges] = histcounts(log10(unitISI)); xlim([0 1000000])
            histogram(unitISI,10.^edges,'Normalization','probability'); %histogram of ISI distribution
            set(gca,'xscale','log') %set log-scale for x-axis
            xticks([10^0 10^2 10^4])
            title(sprintf('unit #%d - log-ISI',n)); ylabel('Probability')
    end

    linkaxes(spwf,'x'); linkaxes(spisi,'x'); %link x-axes in each column of plots
    subplot(spM,spN,(spM*spN)-1); xlabel('Time (ms)');
    subplot(spM,spN,(spM*spN)); xlabel('Time (ms)');

end
