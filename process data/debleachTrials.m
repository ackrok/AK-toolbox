function [out_dff, out_z, time] = debleachTrials (signal, Fs, trialWin)
% Iterative trial-based debleaching photometry analysis
%
% [out_dff, out_z, time] = debleachTrials (signal, Fs, trialWin);
%
% Description: Iterative trial-based debleaching photometry analysis.
% Will calculate baseline in short moving windows, then center
% and normalize these windows, and repeat calculations for temporally
% offset windows within same session (n = #trials). 
% Z-score calculated using mean and stdev of signal spanning entire session
% to minimize confounds in behavior contaminating trial-based analysis.
%
% INPUTS
%   'signal' - processed photometry signal (column vector)
%   'Fs'     - sampling frequency in Hz (scalar)
%   'trialWin' - trial start and end times in samples (matrix, nTrial x 2)
%
% OUTPUTS
%   'out_dff' - output matrix with photometry dF/F by trial
%   'out_z'   - output matrix with photometry z-score by trial
%   'time'    - time vector in seconds relative to trial onset
% 
% Written by Anya Krok, March 2026
%
%% Parameters
blockWindow = 0.5; % moving window duration, in seconds
blockOffset  = 0; % offset, in seconds

trialOn = trialWin(:,1); % trial onset, in samples
if isnan(trialOn(1)); trialOn(1) = 1; end 
trialOff = trialWin(:,2); % trial onset, in samples
maxTrialLen = max(trialOff - trialOn)+1; % maximum trial duration, in samples
nTrial = length(trialOn); % number of trials

%% Correct error values
signal(1) = signal(2); % clip values in first second of session

%% Session-wide z-score parameters
sig_mu = mean(signal, 'omitnan');
sig_std = std(signal, 0, 'omitnan');
if sig_std == 0; sig_std = 1; end % avoid dividing by zero

%% dF/F for entire session
% use baseline FP to identify baseline of specified window using a
% moving window and finding values within a specified percentile
interpType = 'linear'; % 'linear' 'spline' 
fitType = 'interp';  % Fit method 'interp' , 'exp' , 'line'
basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline
%Note: Lower prc are used because the mean of signal is not true baseline
winSize = 10; % Window size for baselining in seconds
winOv = 0; %Window overlap size in seconds
[~, baseline] = baselineFP (signal,interpType,fitType,basePrc,winSize,winOv,Fs);

dFF = (signal - baseline) ./ baseline;

%% Preallocate output: samples x numOffsets x trials
out_dff = nan(maxTrialLen, nTrial); % preallocate output
out_z = out_dff; 
time = makeTime(maxTrialLen, Fs); 
for k = 1:nTrial
    thisOn = trialOn(k);   % onset for this trial
    thisOff = trialOff(k); % offset for this trial
    thisLength = thisOff - thisOn + 1; 
    block = signal(thisOn:thisOff);
    [~,block_base] = baselineFP(block,interpType,fitType,basePrc,blockWindow,blockOffset,Fs);
    block_dff = (block - block_base) ./ block_base;
    block_z = (block_dff - sig_mu) ./ sig_std;
    out_dff(1:thisLength, k) = block_dff;
    out_z(1:thisLength, k) = block_z;
end

%%
% fig = figure; theme(fig, 'light');
% subplot(2,1,1);
% shadederrbar(time, mean(out_dff,2,'omitnan'), std(out_dff,0,2,'omitnan')./sqrt(maxTrialLen), 'k');
% xlabel('time from LED on (s)');
% ylabel('photometry (/deltaF/F)');
% subplot(2,1,2);
% shadederrbar(time, mean(out_z,2,'omitnan'), std(out_z,0,2,'omitnan')./sqrt(maxTrialLen), 'k');
% xlabel('time from LED on (s)');
% ylabel('photometry (z-score)');

