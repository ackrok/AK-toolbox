%% Process Parameters for Photometry Analysis
%Created By: Pratik Mistry
%Created On: 31 January 2019
%Edited On: April 2026, Anya Krok
%
%Description: This is a script with different variables whose values can be
%adjusted depending on the photometry signal that is being processed. The
%name of this file can and should be changed depending on the method and
%GECI used. Please read comments associated with variable
%
%
%% General Parameters
params.acqFs = 5000; % Acquisition sampling frequency
params.dsType = 2; % 1 = Bin Summing; 2 = Bin Averaging; 3 = Traditional (NOT RECOMMENDED)
params.dsRate = 100; % Downsampling rate if you want to downsample the signal
%This dsRate will also be applied to all signals during the analysis
%pipeline
%dsRate = 1 if you do not want to downsample

%% Filter Parameters
params.FP.lpCut = 10; % Cut-off frequency for filter
params.FP.filtOrder = 8; % Order of the filter

%% Baseline Parameters
params.FP.basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline points
%Note: Lower percentiles are used because the mean of signal is not true
%baseline
params.FP.winSize = 10; % Window size for baselining in seconds
params.FP.winOv = 0; %Window overlap size in seconds
params.FP.interpType = 'linear'; % 'linear' 'spline' 
params.FP.fitType = 'interp'; % Fit method 'interp' , 'exp' , 'line'

%% Demodulation Parameters
%When demodulating signals, the filter creates edge artifacts. We record
%for a few seconds longer, so we can remove x seconds from the beginning
%and end
%Adjust the variable to "0" if it's a normal photometry recording
params.FP.sigEdge = 60; %Time in seconds of data to be removed from beginning and end of signal
params.FP.modFreq = 319;

%% Behavior Parameters Parameters


%% Peak Analysis


%% Cross-Correlations


%% Opto-Pulse Analysis
