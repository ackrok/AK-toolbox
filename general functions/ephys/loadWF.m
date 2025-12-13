function WF = loadWF(data,datPath,varargin)
%LoadWF - Load waveforms
%
%   WF = loadWF(data,datPath);
%   WF = loadWF(data,datPath,params);
%
%   Description: This code utilizes organized cluster data (with max channel
%   information) to pull the mean waveforms of all the spikes from the raw
%   data and it stores a randomly selected number of user specified
%   waveforms.
%
%   Input:
%   - data - Data structure containing cluster information
%   - datPath - A string containing the path to the dat file
%   - params (optional) - A structure containing parameters necessary for analysis
%
%   Output:
%   - WF - structure containing waveform information
%
% Author: Pratik Mistry, Tritsch Lab, 2020
% Updated: Anya Krok, 2020

%
% Pull General Information
Fs = data.gen.acqFs;
datName = data.gen.datName;
nChan = data.gen.nChan;
clusters = data.clusters;

% Pull WF params
switch nargin
    case 3
        params = varargin{1}; wfWin = params.WF.window; nWF = params.WF.nToRead;
    case 2
        wfWin = [-3.5 3.5]; % Default: read waveform 3s before and after spike 
        nWF = 100; % Default: read 100 random waveforms 
end

%%
wfWin = ceil((wfWin/1000)*Fs);
wfTime = wfWin(1):wfWin(2); wfTime = (wfTime'/Fs)*1000;
wfDur = length(wfTime);

% Create Memory Map to Dat File
datFF = fullfile(datPath,datName);
datStruct = dir(datFF);
nBytes = numel(typecast(cast(0,'int16'),'uint8'));
nSamp = datStruct.bytes/(nChan*nBytes);
dataMap = memmapfile(datFF, 'Format', {'int16', [nChan nSamp], 'dat'});

nClu = length(clusters);
WF(1:nClu) = struct('top',zeros(wfDur,nWF),'topMu',zeros(wfDur,1),'topSem',zeros(wfDur,1),'time',zeros(wfDur,1));

%The following for loop will go through all the clusters and pull the
%waveforms from the specified spike times within a userset window. It also
%creates a vector containing random integers selected from 1 -> # of spikes
%to pull the individual waveforms.
tic
h = waitbar(0,'Loading Waveforms from DAT File');
for n = 1:nClu
    clustID = clusters(n).clusterID; % Pull the cluster ID
    tmpST = clusters(n).spikeTimes; % Pull the spike times from that cluster
    tmpST = int64(tmpST * Fs);
    nSpikes = length(tmpST); % Store number of spikes
    tmpWF = zeros(nSpikes,wfDur); % Initialize a matrix that will house all the waveforms from a cluster
    randInd = randi(nSpikes,nWF,1); % Create a vector of randomly selected integers from 1 -> # Spikes
    maxChan = clusters(n).maxChannel; % Pull the max channel used when pulling waveforms from raw data
    for m = 1:nSpikes
        ind1 = tmpST(m)+wfWin(1); %Pull index 1 from spike time
        ind2 = tmpST(m)+wfWin(2); %Pull index 2 from spike time
        %The following if statement checks to see if index 1 is below 1 and
        %if index 2 is past the total number of samples in a dataset
        if ind1 < 1
            ind1 = 1; ind2 = wfDur;
        elseif ind2 > nSamp
            ind1 = nSamp - wfDur + 1; ind2 = nSamp;
        end
        tmpWF(m,:) = double(dataMap.Data.dat(maxChan,[ind1:ind2])); %Store waveforms in temp matrix
    end
    WF(n).top = tmpWF(randInd,:); %Extract random waveform
    WF(n).topMu = mean(tmpWF,1); WF(n).topMu = WF(n).topMu(:); %Extract mean waveform from all spikes
    WF(n).topSem = SEM(tmpWF,1); WF(n).topSem = WF(n).topSem(:); %Extract standard error of mean from all spikes
    WF(n).time = wfTime(:);
    waitbar(n/nClu,h);
    %fprintf('Loaded Waveform for Cluster %d.... Completed %d out of %d\n',clustID,n,nClu);
end
disp('Finished Loading Waveforms!'); close(h);
toc
end
