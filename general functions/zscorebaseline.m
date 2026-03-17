function zscoresignal = zscorebaseline(t, signal, baseline_sec)
% zscorebaseline
%   Z-score a signal relative to a baseline window.
%   If baseline mean is negative or very small, shift the signal so that
%   the baseline mean becomes a positive reference before z-scoring.
%
% INPUTS:
%   t           : time vector (s), same length as signal
%   signal      : vector (e.g., dFF)
%   baseline_sec: scalar, length of baseline window in seconds (e.g. 600)
%
% OUTPUT:
%   zscoresignal: z-scored signal

% ---- define baseline mask (0..baseline_sec) ----
baseline_mask = t >= 0 & t <= baseline_sec;

if ~any(baseline_mask)
    error('zscorebaseline: no points found in baseline window 0–%.1f s.', baseline_sec);
end

% ---- compute baseline mean ----
base_vals    = signal(baseline_mask);
mu_base_raw  = mean(base_vals);

% ---- optional: adjust signal if baseline is negative or too small ----
% This is analogous to aligning a model's mean to the raw baseline mean.
% Here we want a baseline mean that is clearly positive to avoid
% interpretation issues when later using dFF or for plotting.
if mu_base_raw <= 0 || abs(mu_base_raw) < 1e-6
    % Choose a simple positive offset based on the raw baseline magnitude
    offset        = abs(mu_base_raw) + 1e-6;
    signal_shift  = signal + offset;
    base_vals     = signal_shift(baseline_mask);
    mu_base       = mean(base_vals);
else
    signal_shift = signal;
    mu_base      = mu_base_raw;
end

% ---- compute baseline std ----
sd_base = std(base_vals);

if sd_base <= 0 || isnan(sd_base)
    error('zscorebaseline: baseline standard deviation is zero or NaN.');
end

% ---- z-score using adjusted signal and baseline stats ----
zscoresignal = (signal_shift - mu_base) / sd_base;

end