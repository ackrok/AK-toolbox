function eventTimes = adjEvtTimes(eventTimes,tBefore,tAfter,totSamp)
evtAfter = eventTimes + tAfter; evtBefore = eventTimes + tBefore;
logInd = evtAfter < totSamp | evtBefore > 1;
eventTimes = eventTimes(logInd);
end