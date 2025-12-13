function [spk] = getSpikeProp(clusters)
%Cluster temporal properties from spike times
% 
% Description: This function will use cluster spike times to compute 
%   temporal/spiking properties for each cluster.
%
% [spk] = getSpikeProp(clusters)
%
% INPUT
%   'clusters' - structure containing cluster spike time information
%
% OUTPUT
%   'spk.fr'    - (Hz) mean firing rate 
%   'spk.ISI'   - (*ms*) difference between successive spikeTimes
%   'spk.pISI2' - proportion of time associated with long ISI (> 2 sec)
%   'spk.CV'    - coefficient of variance: CV = 0 regular, CV = 1 random, CV > 1 bursting
%   'spk.burst' - structure with burst properties
%   'spk.bprop' - proportion of spikes within bursts
%
% Created By: Anya Krok, June 2019

%% Initialize Variables
fr          = {};
ISI         = {}; 
pISI2       = {}; 
CV          = {}; 
burst       = {};
bprop       = {};

stAll = {clusters.spikeTimes}; %Extract all spike times into cell array
lastSt = cellfun(@(x) x(end), stAll); %Vector of last spike for each unit, in seconds
stLast = max(lastSt); %Time of last spike over all units, in seconds 

%% Temporal Properties
for n = 1:length(clusters)
% Extract spike times from data structure
    st = double(clusters(n).spikeTimes);
    
% Firing rate
    fr{n} = 1/mean(diff(st));
    %fr{unitNum} = length(data.clusters(unitNum).spikeTimes)/data.final.time(end);
    
% ISI: difference between successive spike times
    tmpISI = diff(st); 
    ISI{n} = tmpISI.*1000; %convert to ms
    
% pISI2: determine proportion of time associated with long ISIs (ISI > 2s),
    % used to distinguish between phasically and tonically active neurons
    pISI2{n} = sum(tmpISI(tmpISI > 2));     %Sum ISIs longer than 2s 
    pISI2{n} = pISI2{n}/stLast;             %then divide sum by total time
    
% Coefficient of Variance, CV = 0 regular, CV = 1 random, CV > 1 bursting
    CV{n} = std(tmpISI)/mean(tmpISI); 

% Bursting: total number, duration, mean IBI, #spikes/burst, onset, offset
    tmpBurst = getBurst(st);
    burst{n} = tmpBurst;
    bprop{n} = tmpBurst.prop; % Extract proportion of spikes in burst
    %unless otherwise specified, bursting is 3+ spikes with ISI <50ms
end

%% Generate Output Structure
spk = struct('fr',fr,'ISI',ISI,'pISI2',pISI2,'CV',CV,'burst',burst,'bprop',bprop);

end