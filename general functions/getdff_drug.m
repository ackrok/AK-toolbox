function [dff, z] = getdff_drug(signal, trend, Fs, win)
% getdff_drug  
% Compute deltaF/F and z-score for photometry signal after detrending
%
% Syntax:
%   [dff, z] = getdff_drug(signal, trend, Fs, win)
%
% Inputs:
%   signal  - photometry signal prior to detrending
%   trend   - fitted trend of data, from detrendFP_drug fxn
%   Fs      - sampling rate, in Hz (scalar)
%   win     - time point, in seconds, that defines end of baseline
%               will be used to compute mask = 1 : win*Fs
%
% Outputs:
%   dff - deltaF/F: (detrend - F0) ./ F0  (vector)
%   z   - z-score:   (detrend - F0) ./ sigma (vector)
%
% Written by Anya Krok, March 2026

% Detrend photometry signal
detrend = signal(:) - trend(:);

% Compute baseline sample indices: 1 through win*Fs
nBaseline = round(win * Fs);      % number of baseline samples
nBaseline = max(1, min(nBaseline, numel(detrend))); % clamp to valid range
mask = 1 : nBaseline;             % baseline mask (indices)

% Baseline mean (F0) computed from detrended baseline region
F0 = mean(detrend(mask));

% If baseline mean is negligibly small, add DC offset estimated from raw signal
% This avoids dividing by near-zero and preserves relative scale.
if F0 <= 1e-3
    offset  = mean(signal(mask)); % DC offset
    detrend = detrend + offset;
    F0 = F0 + offset;
end
% Baseline standard deviation used for z-scoring
sigma = std(detrend(mask));

% Compute deltaF/F (fractional change relative to baseline mean)
% Note: if F0 == 0 this will produce Inf/NaN; upstream code should check or avoid zero baseline.
dff = (detrend - F0) ./ F0;

% Compute z-score relative to baseline mean and baseline standard deviation
% Note: if sigma == 0 this will produce Inf/NaN; consider handling this case if needed.
z = (detrend - F0) ./ sigma;