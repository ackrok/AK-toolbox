function [burst] = getBurst(st, varargin)
%Burst characterization
%
% Description: This function will use spike times for a given unit to
%   characterize burst events, defined as 3+ spikes with ISI < 50ms
%
% [burst] = getBurst(st, varagin)
%
% INPUT
%   'st' - vector of spike times for a single unit, in seconds
%   varargin
%       {1} 'threshold' - threshold ISI value, default 0.05 seconds
%       {2} 'thresSpikes' - (#) number of spikes, default 3+ spikes
%
% OUTPUT
%   'burst' - structure with sub-structures:
%       'burst.num'     - total number of bursts
%       'burst.prop'    - proportion of spikes that are within bursts
%       'burst.muDur'   - (*ms*) average burst duration
%       'burst.muIBI'   - (*ms*) average inter-burst interval
%       'burst.muSpk'   - average number of spikes in a burst
%       'burst.on'      - (*sec*) burst onset 
%       'burst.off'     - (*sec*) burst offset
%
% Created By: Anya Krok, June 2019

%% Variables
ISI = diff(st); % ISI: difference between successive spike times

if nargin == 1
    threshold = 0.050;   %threshold value to define bursting, in seconds
    thresSpk = 3; %threshold for number of spikes necessary for burst
        %note that >3 spikes with ISI < 50 corresponds to >2 intervals between
        %spikes so we use thresSpikes-1 for finding regions of 3+ spikes
elseif nargin == 3
    threshold = varargin{1};
    thresSpk = varargin{2};
end

%%
% Find logical vector where ISI < threshold    
binaryVector = ISI < threshold; 

% Label each region with a label/identifier number
[labelVector, numRegions] = bwlabel(binaryVector);

% Measure lengths of each region and the indexes
measurements = regionprops (labelVector, ISI, 'Area', 'PixelValues', 'PixelIdxList');

% Find regions where 3+ spikes and put values into a cell of cell array
burstArray = cell(1,numRegions); burstIdx = cell(1,numRegions);
for k = 1:numRegions
    if measurements(k).Area >= thresSpk
        burstArray{k}   = measurements(k).PixelValues;
        burstIdx{k}     = measurements(k).PixelIdxList; %idx in spikeTimes
    else
        burstArray{k} = []; burstIdx{k} = [];
    end
end

% Remove empty cells from cell arrays
burstArray  = burstArray(~cellfun('isempty',burstArray));
burstIdx    = burstIdx(~cellfun('isempty',burstIdx)); 
burstNum    = length(burstArray);

% Characterize bursts: burst duration, spike number per burst, IBI
allBurstDur = []; allBurstSpikes = []; allIBI = []; 
burstOn = []; burstOff = [];

for n = 1:burstNum
    %sum of ISI from first index of current burst to last index of burst
    allBurstDur(n) = sum(burstArray{1,n}); 
    
    %number of spikes = 1 + number of intervals between spikes
    allBurstSpikes(n) = 1 + length(burstArray{1,n});
    
    %IBI is sum of inter-spike intervals between first index of current
    %burst and last index of previous burst, excluding those indices
    if n > 1; prevBurst = burstIdx{1,n-1}; 
        else; prevBurst = 0;
    end
    thisBurst = burstIdx{1,n};
    allIBI(n) = sum(ISI(prevBurst(end)+1 : thisBurst(1)-1));
    
    %burst onset and offset as index: relate to time as spikeTimes(burstOn)
    burstOn(n) = thisBurst(1);
    burstOff(n) = thisBurst(end);
end

burst.num   = burstNum;
burst.prop  = sum(allBurstSpikes)/length(st); 
burst.muDur = mean(allBurstDur); 
burst.muSpk = mean(allBurstSpikes); 
burst.muIBI = mean(allIBI);      
%burst.semDur     = SEM(allBurstDur, 2);
%burst.semSpikes  = SEM(allBurstSpikes, 2);
%burst.semIBI     = SEM(allIBI, 2);
burst.on    = st(burstOn); 
burst.off   = st(burstOff);

end