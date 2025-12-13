function splitDat(totalChan,timeWin,chan2pull,Fs,memToUse)
%splitDat - Split Dat Files
%
%   Usage:
%       splitDat(totalChan,timeWin,chan2pull,Fs,memToUse)
%       splitDat(totalChan,{},chan2pull,Fs,memToUse) - If splitting channels over entire recording
%       splitDat(totalChan,{t,'end'},chan2pull,Fs,memToUse) - If splitting
%           channels from time point t seconds to end of recording
%       splitDat(totalChan,{'start',t2},chan2pull,Fs,memToUse) - If splitting
%           channels from the beginning to a time point t2
%       splitDat(totalChan,[t1,t2],chan2pull,Fs,memToUse) - If splitting
%           channels from time point t1 to time point t2, in seconds
%
%   Description:
%       This function will parse .dat file based on specified channels to
%       extract and/or beginning and end time points.
%
%   Input:
%       totalChan - Total number of channels in dat file
%           * Neuronexus Buzsaki A32 probe has 43 total channels
%           * Masmanidis 128ch probe has 142 total channels
%       timeWin - A vector or cell array containing the time in seconds to pull from
%       or the values 'start' and 'end' to indicate the beginning or end
%       of a recording
%       chan2pull - A vector containing channels to extract
%           * Neuronexus Buzsaki A32 probe has 3 aux channels
%                   Probe channels are [1:32], ADC channels are [36:43]
%           * Masmanidis 128ch probe has 6 aux channels
%                   Probe channels are [1:128], ADC channels are [135:142]
%       Fs - Acquisition Sampling Rate
%       memToUse - Amount of memory in GB to use while parsing the file
%
%   Output:
%       The function writes a new file to the same path as the file being
%       parsed
%
%
%   Author: Pratik Mistry, Anya Krok, Tritsch Lab, 2020

[datName,datPath] = uigetfile('*.dat','Select File to Parse'); %Get the original dat file 
datFF = fullfile(datPath,datName); %Create file full path
while(1)
    datInput = inputdlg('Enter Name for New File'); %Ask the user to input new name for file
    datName_new = datInput{1};
    [~,datNewExt] = strtok(datName_new,'.');
    if (isempty(datNewExt))
        datName_new = [datName_new,'.dat'];
    end
    if(~strcmp(datName_new,datName))
        break;
    else
        err = msgbox('Please create a filename different than the original');
        uiwait(err);
    end
end
try
    datMap = getDatMemMap(datFF,totalChan);  %Create a memory map to the dat file -- allows for faster access
catch
    msgbox('Incorrect number of channels!');
    return;
end
totSamp = datMap.Format{2}(2); %Pull the total number of samples in time from file
[tBegin,tEnd] = parseTimeWin(timeWin,Fs,totSamp);
sampToPull = length(tBegin:tEnd);
buff = (memToUse*(10^9))/(16*length(chan2pull));
nChunks = ceil((sampToPull/buff));

fidout = fopen(fullfile(datPath,datName_new),'w');
fopen(fidout);
for ii = 1:nChunks
    fprintf('Splitting dat file... writing chunk %d out of %d...\n',ii,nChunks);
    ind1 = int32(((ii-1)*buff)+tBegin);
    if ii == nChunks
        ind2 = sampToPull;
    else
        ind2 = int32((ii*buff)+tBegin);
    end
    tmpData = datMap.Data.dat(chan2pull,[ind1:ind2]);
    fwrite(fidout,tmpData,'int16');
end
fclose(fidout);
fprintf('Finished splitting the file!\n');
end

function [tBegin,tEnd] = parseTimeWin(timeWin,Fs,totSamp)
if(isempty(timeWin))
    tBegin = 1; tEnd = totSamp;
elseif iscell(timeWin)
    if ischar(timeWin{1})
        tBegin = 1;
    else
        tBegin = timeWin{1}*Fs;
    end
    if ischar(timeWin{2})
        tEnd = totSamp;
    else
        tEnd = timeWin{2}*Fs;
    end
elseif isvector(timeWin)
    tBegin = timeWin(1)*Fs;
    tEnd = timeWin(2)*Fs;
end
if (tEnd>totSamp)
    tEnd = totSamp;
end
if (tBegin<1)
    tBegin = 1;
end
end