function [muDeltaT, range] = getMuDeltaT(st, events, varargin)
% Computation of mean delta-t of unit cell-pair spikes.
%
% [muDeltaT, range] = getMuDeltaT(st, events, varargin)
% [muDeltaT, range] = getMuDeltaT(st, events, range)
%
% Description: 
%   Identify spike of n'th unit that is closest to specified event time,
%   then find spike of m'th unit closest that of the n'th unit. 
%   Delta-t is the time difference between spike(m) and spike(n).
%   Repeat for all cell pairs, then repeat over multiple event times. 
%   Average delta-t for each cell pair over all events, and then average
%   all the cell pairs for each latency/shift.
%
% INPUT
%   'st' - cell array with unit spike times, in seconds
%   'events' - vector with event times, in seconds
%   'range' (optional) - vector of event latencies, in seconds
%
% OUTPUT
%   'muDeltaT' - matrix with delta-t averaged over all events for each cell
%       pair at each latency, where columns correspond to cell pairs and
%       rows correspond to different latencies from event times.
%   'range' - vector of event latencties, in seconds
%
% Anya Krok, March 2020
% 

%% Input Variables
    units = 1:length(st); nUnits = length(st); 
    switch nargin
        case 3
            range = varargin{1}; 
        case 2
            range = [-5:0.5:5]; %Time vector 
    end

%% deltaT of unit spike time to unit spike time 
    muDeltaT = zeros(length(range), nchoosek(nUnits,2)); %Initialize output matrix
    h1 = waitbar(0, 'Computing delta t for Multiple Event Latencies');
    for a = 1:length(range) %Repeat calculation of delta t for multiple shifts of event time
        onTime = events + range(a); %Shift event time by some number of seconds
        dt_shift = zeros(length(onTime), nchoosek(nUnits,2)); %Initialize output matrix
        for z = 1:length(onTime) %Repeat calculation of delta t for all events
            event = onTime(z); %Event time
            deltaT = nan(nUnits, nUnits); %Initialize output matrix
            for x = 1:nUnits; n = units(x);
                st_n = st{n}; %Extract spike times for n'th unit
                [~, idx_n] = min(abs(st_n - event)); %Find closest spike to event by minimizing difference between st and event
                time_n = st_n(idx_n); 
                for y = 1:nUnits; m = units(y); 
                    st_m = st{m}; %Extract spike times for m'th unit
                    [~, idx_m] = min(abs(st_m - time_n)); %Find closest spike
                    time_m = st_m(idx_m); 

                    deltaT(x,y) = abs(time_m - time_n).*1000; %Absolute time between spikes of two units, in milliseconds
                end
            end
            ii = ones(size(deltaT));
            rmvidx = tril(ii,-1);  %Remove top half with diagonal
            deltaT(~rmvidx) = [];  %Replace deleted elements with empty vector
            dt_shift(z,:) = deltaT; 
        end
        muDeltaT(a,:) = mean(dt_shift,1); %Average across all bouts for each cell pair
        waitbar(a/length(range), h1); 
    end
    close(h1); 
    fprintf('\n Analysis Complete. Computed mean delta-t for unit cell pairs. \n');

end
