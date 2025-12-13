function eventMat = alignTrace2Event(trace,eventTimes,tBefore,tAfter,nSamp)
%alignTrace2Event - Align time series traces to event times
%
%   Usage:
%       eventMat = alignTrace2Event(trace,eventTimes,tBefore,tAfter,nSamp)
%
%   Description: This function creates a matrix of time series data
%   aligned to event times specifed by the user within a certain window
%   
%   Input:
%       trace - Time series trace to be used
%       eventTimes - A vector containing event times in samples
%       tBefore - Amount of time before the event time to capture
%       tAfter - Amount of time after the event time to capture
%       nSamp - Number of samples of the window
%
%   Output:
%       eventMat - Matrix containing event aligned traces
%
%   Author: Pratik Mistry, 2020
%
nEvents = length(eventTimes);
eventMat = zeros(nSamp,nEvents);
for n = 1:nEvents
    tmpEvent = eventTimes(n);
    ind1 = tmpEvent + tBefore; ind2 = tmpEvent + tAfter;
    if ind1<1 || ind2 > length(trace)
        eventMat(:,n) = nan(size(eventMat(:,n)));
        eventTimes(n) = NaN;
    else
        eventMat(:,n) = trace(ind1:ind2);
    end
end
if sum(any(isnan(eventMat),1)) > 0
    eventMat(:,(any(isnan(eventMat),1))) = [];
end
if sum(any(isnan(eventTimes))) > 0; eventTimes(any(isnan(eventTimes))) = []; end
end