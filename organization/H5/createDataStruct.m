function [data] = createDataStruct(wsData,animalName,expDate)
%%Create Data Structure from WaveSurfer to MAT function
%
%   [data] = createDataStruct(wsData,animalName,expDate)
%
%   Description: This file takes a data structure generated from the
%   extractH5_WS function and it organizes the acquired traces into a lab
%   specific format by asking the user specify which type of recording the
%   trace is: (Photometry, Reference, Wheel, Opto Pulse)
%
%   Input:
%   - wsData - Data structure generated from extractH5_WS
%   - animalName - Name of animal/exp to add to data structure
%   - expDate - Date that the experiment took place
%
%   Output:
%   - data - New data structure that is cleaner to parse
%
%
%   Author: Pratik Mistry 2019

data = initDS; %Intialize data structure
data.mouse = animalName; data.date = expDate; %Add mouse name and experiment date to structure
nSweeps = length(wsData.sweeps);
%The following for-loop will go through all the sweeps and organize the
%data within each sweep depending on the trace name
traceNames = wsData.allNames;
nTraces = wsData.header.nTraces;
data.gen.acqFs = wsData.header.AcquisitionSampleRate; %Pull Sampling Rate
data.gen.startTime = wsData.header.ClockAtRunStart; %Pull time at start
for x = 1:nTraces
    tName = traceNames{x};
    tName(find(tName==' ')) = [];
    traceNames{x} = tName;
    %Following line asks the user to select the type for the trace
    switch tName % AK added 220921
        case {'ACh','DA','dLight','rACh','rDA','5-HT'}
            choice = 1;
        case {'Ref_ACh','Ref_DA','ref_ACh','ref_DA','ref_green','ref_red','Ref_green','Ref_red'}   
            choice = 2;
        case {'Wheel','wheel'}
            choice = 3;
        case {'Reward','reward'}
            choice = 5;
        case {'Lick','lick'}
            choice = 6;
        otherwise    
            choice = menu(['Select an option for trace: ',tName],'Photometry','Reference',...
                'Wheel Encoder','Opto Pulses','Reward Behavior','Lick','Camera Trigger','Delete');
    end
    traceType(x) = choice;
end
for sweepNum = 1:nSweeps
    Ls = size(wsData.sweeps(sweepNum).allData,1);
    FPind = 1; RefInd = 1; %Need to initialize index for FP and Ref sigs because there may be multiple
    RewInd = 1; LickInd = 1; %Initialize index for Reward, Lick sigs because may be multiple
    for n = 1:nTraces
        %Following line asks the user to select the type for the trace
        choice = traceType(n);
        %The switch statement will organize the traces into specific fields
        %in the new data structure depending on user input
        switch choice
            case 1
                data.acq(sweepNum).FP{FPind,1} = wsData.sweeps(sweepNum).allData(:,n);
                data.acq(sweepNum).FPnames{FPind,1} = traceNames{n};
                FPind = FPind+1;
            case 2
                data.acq(sweepNum).refSig{RefInd,1} = wsData.sweeps(sweepNum).allData(:,n);
                data.acq(sweepNum).refSigNames{RefInd,1} = traceNames{n};
                RefInd = RefInd+1;
            case 3
                data.acq(sweepNum).wheel = wsData.sweeps(sweepNum).allData(:,n);
            case 4
                data.acq(sweepNum).opto  = wsData.sweeps(sweepNum).allData(:,n);
            case 5
                data.acq(sweepNum).rew{RewInd,1} = wsData.sweeps(sweepNum).allData(:,n);
                data.acq(sweepNum).rewNames{RewInd,1} = traceNames{n};
                RewInd = RewInd+1;
                % data.acq(sweepNum).rew   = wsData.sweeps(sweepNum).allData(:,n);
            case 6
                data.acq(sweepNum).lick{LickInd,1} = wsData.sweeps(sweepNum).allData(:,n);
                data.acq(sweepNum).lickNames{LickInd,1} = traceNames{n};
                LickInd = LickInd+1;
                % data.acq(sweepNum).lick  = wsData.sweeps(sweepNum).allData(:,n);
            case 7
                data.acq(sweepNum).cam   = wsData.sweeps(sweepNum).allData(:,n);
            case 8
                
        end
    end
    data.acq(sweepNum).nFPchan = length(data.acq(sweepNum).FP);
    timeVec = [1:Ls]/data.gen.acqFs;
    data.acq(sweepNum).time = timeVec';
end
end

function data = initDS()
data = struct('mouse',[],'date',[],'acq',struct());
data.acq = struct('FPnames','','nFPchan',[],'FP',[]);
end