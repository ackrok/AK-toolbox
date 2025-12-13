function [pause] = getPause(st, varargin)
%Pause characterization
%
% Description: This function will use spike times for a given unit to
%   characterize pauses
%
% [pause] = getPause(st, varargin)
%
% INPUT
%   'st' - vector of spike times for a single unit, in seconds
%   varargin
%       {1} 'coreInt' - length of interval for identifying pauses (*sec*)
%       {2} 'minPause' - minimal length of final pause (*sec*)
%       {3} 'addLimit' - upper limit of added intervals
%
% OUTPUT
%   'pause' - structure with sub-structures:
%       'pause.num'     - total number of pauses
%       'pause.muDur'   - (*ms*) average pause duration
%       'pause.muIPI'   - (*ms*) average inter-pause interval
%       'pause.on'      - (*sec*) pause onset
%       'pause.off'     - (*sec*) pause offset
%
% Created By: Anya Krok, June 2019

%% Input Variables
ISI = diff(st); % ISI: difference between successive spike times

%coreInt = 0.25; %Length of interval for identifying pause, in seconds
coreInt = 10.*median(ISI); %Can be determined according to firing rate of cell

addLimit = 5; %Upper limit of added intervals, from each side of the core interval
addRange = [-1.*addLimit : addLimit];

minPause = 0.3; %minimal length of the final pause, in seconds
%minPause can also be mean onset of second peak in ISI histogram of pausers
%who have bimodal ISI histograms

spikeMerge = 1; %maximal number of sikes enabling merging of pauses
%value of X indicates that two pauses separated by X spikes or less can be
%merged 

params = struct;
params.coreInt = coreInt; 
params.minPause = minPause;
params.addLimit = addLimit; 

%% Output Variables
pause = struct;

%% Pause Detection
detectIdx = find(ISI > coreInt); %Detect ISIs longer than a certain length (defined by core interval)
if isempty(detectIdx)
    return
end

fr = 1/mean(ISI); %Average firing rate of spike train, in 1/seconds
PSfun = @(x) fr*x*exp(fr*x);
%Poisson surprise value (S) of the pause: S = -logP(n), where P(n) is
%probability of finding n spikes or less in an interval of T seconds
x = detectIdx(1,1);
for n = 2:length(detectIdx)-1
    if detectIdx(n+1,1) > detectIdx(n,1) + addLimit
        x = [x, detectIdx(n+1)];
    end
end
detectIdx = x; %Indices of ISI intervals satisfying pause criteria


%% Pause Detection
detectIdx = find(ISI > coreInt); %Find indices where ISI > threshold
if isempty(detectIdx)
    return
end

fr = 1/mean(ISI); %Average firing rate of spike train, in 1/seconds
PSfun = @(x) fr*x*exp(fr*x); %Poisson surprise value (S) of the pause
% S = -logP(n), where P(n) is probability of finding n spikes or less in an interval of T seconds

pauseAll = struct; 
for x = 1:length(detectIdx) 
    idx = detectIdx(x); %Index of ISI
    T   = ISI(idx);       %Length of original ISI, in seconds
    S   = PSfun(T);       %Poisson surprise value
   
    pauseAll(x).idx = idx;   
    pauseAll(x).dur = sum(ISI(pauseAll(x).idx)); 
    
    for y = addRange %We will test intervals before and after core interval
        idxTest = idx + y; 
        if idxTest <= length(ISI) && idxTest >= 1 && idxTest ~= idx 
            Tnew = T + ISI(idxTest); %New period is original plus added ISI 
            Snew = PSfun(Tnew);      %New Poisson surprise value for new period
            if Snew > S     %IF Poisson surprise value increases
               T = Tnew;    %THEN add additional interval to core interval
               pauseAll(x).idx = sort([pauseAll(x).idx, idxTest]); 
               pauseAll(x).dur = sum(ISI(pauseAll(x).idx)); 
            end
        end
    end
end

%% Combine overlapping pauses    
for x = 1:length(pauseAll)
    for y = 1:length(pauseAll)
        check = any(ismember(pauseAll(x).idx, pauseAll(y).idx)); %Determine if there is any overlap in indices
        if check
            idx = sort(unique([pauseAll(x).idx, pauseAll(y).idx])); %Combine indices
            pauseAll(x).idx = idx;
            pauseAll(x).dur = sum(ISI(pauseAll(x).idx)); 
        end
    end
end

%% 
uniPause = unique([pauseAll.dur]); %Identify unique pause duration times
uniPause = uniPause(uniPause > minPause); %Minimal length of the final pause

pauseUni = struct;
for x = 1:length(uniPause)
    y = find(ismember([pauseAll.dur],uniPause(x)),1); %Index of first match to unique pause duration time
    pauseUni(x).idx = pauseAll(y).idx; %Save unique index into new structure
    pauseUni(x).dur = sum(ISI(pauseUni(x).idx)); 
    pauseUni(x).on  = st(pauseUni(x).idx(1)); %Start time of pause, in seconds
    pauseUni(x).off = st(pauseUni(x).idx(end)+1); %End time of pause in seconds
end
tab = struct2table(pauseUni); %Covert struct array to a table
sortTab = sortrows(tab,'on'); %Sort table by 'on'
pauseUni = table2struct(sortTab); %Change back into struct array
for x = 2:length(pauseUni)
    pauseUni(x).IPI = pauseUni(x).on - pauseUni(x-1).off; %Inter-pause interval, in seconds
end

%% Pause properties
pause.num = length(pauseUni);       %Number of pauses in recording
pause.muDur = mean([pauseUni.dur]); %Average duration of pause, in seconds
pause.muIPI = mean([pauseUni.IPI]); %Average inter-pause interval, in seconds
pause.on = [pauseUni.on];        %Start time of pauses, in seconds
pause.off = [pauseUni.off];      %End time of pauses, in seconds
pause.params = params;

end