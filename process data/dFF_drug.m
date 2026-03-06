function out = dFF_drug(signal, Fs, winBase)
%dFF_drug - Baseline adjust photometry signal to get dF/F, utilizing
%pre-drug baseline period
%
%   [out] = dFF_drug(signal, Fs, baseWin)
%
%   Description: This code will baseline adjust the photometry signal using
%   baseline period, prior to drug administration
%
%   Input:
%   - signal - Photometry signal to baseline
%   - winBase - Window for baseline period, in seconds (eg, [0 600])
%   - Fs - Sampling Rate, in Hz
%
%   Output:
%   - out
%
%   Author: Anya Krok, March 2026
%   adapted from baseline FP
%

%Ensure the FP vector is a column vector instead of row vector --> Faster
%computation
signal = signal(:);
winBase = winBase * Fs; %Convert window from seconds to samples
baseline = signal( winBase(1) : winBase(2) ); 

Ls = length(signal); %Get length of photometry trace
L = 1:Ls; L = L'; %Create a column vector from 1 to total data points in trace
nPts = floor(Ls/(winSize-winOv)); %Determine number of baseline points to be found

%Baseline 
% (1) use baseline FP to identify baseline of specified window using a
%     moving window and finding values within a specified percentile
interpType = 'linear'; % 'linear' 'spline' 
fitType = 'interp';  % Fit method 'interp' , 'exp' , 'line'
basePrc = 50; % Percentile value from 1 - 100 to use when finding baseline
%Note: Lower prc are used because the mean of signal is not true baseline
winSize = 10; % Window size for baselining in seconds
winOv = 0; %Window overlap size in seconds
[~, outbaseline, ~] = baselineFP (baseline,interpType,fitType,basePrc,winSize,winOv,Fs);
meanbaseline = mean(outbaseline);
% (2) mean of baseline
meanbaseline = mean(baseline);

% Computer dF/F of entire signal
dFF = (signal - meanbaseline) ./ meanbaseline;
dFF = dFF * 100; 

end
