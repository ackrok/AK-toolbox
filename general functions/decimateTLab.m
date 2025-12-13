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