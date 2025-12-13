function [lfp_full] = loadLFP(nChan,chan2pull,acqFs,dsFs,varargin)
%loadLFP - Load LFPs from In Vivo Ephys Recordings
%
%   Usage:
%       lfp = loadLFP(nChan,chan2pull,acqFs,dsFs)
%
%   Description: This function pulls a single channel recording from an In
%   Vivo Ephys recording and decimates the trace to allow for LFP analysis
%   
%   Input:
%       nChan - Number of channels in the dat file
%       chan2pull - The single channel to pull from a recording
%       acqFs - The raw acquisition sampling frequency
%       dsFs - The new downsampled frequency
%       datFF - fullfile name for dat file to extract LFP from
%
%   Output:
%       lfp - The decimated channel data
%
%
%   Author: Pratik Mistry, 2020
%

switch nargin
    case 4
        [datFile,datPath] = uigetfile('*.dat','Select Dat File to Extract Single Channel LFP');
        datFF = fullfile(datPath,datFile);
    case 5
        datFF = varargin{1};
end
    
datMap = getDatMemMap(datFF,nChan);
lfp_full = double(datMap.Data.dat(chan2pull,:));
% lfp = decimateTLab(lfp_full,acqFs,acqFs/dsFs);
end


function data = decimateTLab(rawData,rawFs,dsRate)
%decimate - Reduce sampling rate of data
%
%   [data] = decimateTLab(rawData,rawFs,dsRate)
%
%   Description: This function reduces the sampling rate of a data set by
%   applying a 10-th order butterworth filter at a cut-off frequency that
%   is 45% of the new sampling rate. And then downsamples the data by
%   picking every nth data point.
%
%   Input:
%   - rawData - Original raw data set
%   - rawFs - Original sampling rate
%   - dsRate - Downsampling rate
%
%   Output:
%   - data - Decimated data set
%
%   Author: Pratik Mistry, 2019
%
    Fs = rawFs/dsRate;
    lpFilt = designfilt('lowpassiir', 'FilterOrder', 10, 'HalfPowerFrequency', floor(Fs*0.45), 'SampleRate', rawFs);
    data = filtfilt(lpFilt,rawData);
    data = data(1:dsRate:end);
end