function [corrected_signal, dFF, metrics] = iso_modeA(signal, fs)
% iso_modeA
%   Mode A: Use a 3-point quadratic bleaching model.
%   - Clamp first 2 s
%   - Take 0–600 s baseline from signal
%   - Optionally smooth baseline (currently off)
%   - Build a quadratic that passes through (t_start, t_mid, t_end)
%   - Evaluate model over full recording
%   - Subtract model from signal
%   - Compute ΔF/F using same baseline
%
% INPUTS:
%   signal : [N x 2], col1 = time (s), col2 = signal (e.g., red)
%   fs     : sampling rate (Hz)
%
% OUTPUTS:
%   corrected_signal : [N x 1] corrected signal
%   dFF              : [N x 1] ΔF/F
%   metrics          : struct with QC metrics

t     = signal(:,1);
F_sig = signal(:,2);

% ----- clamp first 2 s -----
clamp_sec = 2;
[~, idx2] = min(abs(t - clamp_sec));
idx_early = t < clamp_sec;
F_sig(idx_early) = F_sig(idx2);

% ----- define baseline window (0–600 s) -----
baseline_sec = 600;
baseline     = t >= 0 & t <= baseline_sec;
if ~any(baseline)
    error('iso_modeA: no samples in baseline window 0–%d s.', baseline_sec);
end
t_base = t(baseline);
F_base = F_sig(baseline);

% ----- lightly smooth baseline to reduce noise (currently no smoothing) -----
smooth_win_sec     = 0;                        % 0 s -> no smoothing
smooth_win_samples = max(1, round(smooth_win_sec * fs));
F_base_smooth      = movmean(F_base, smooth_win_samples);

% ----- build 3-point quadratic bleaching model -----
t_start = t(1);
t_mid   = min(600, t(end));   % ensure mid is within data
t_end   = t(end);

[~, i_start] = min(abs(t - t_start));
[~, i_mid]   = min(abs(t - t_mid));
[~, i_end]   = min(abs(t - t_end));

y_start = F_sig(i_start);
y_mid   = F_sig(i_mid);
y_end   = F_sig(i_end);

t_pts = [t_start; t_mid; t_end];
y_pts = [y_start; y_mid; y_end];

% fit a quadratic that passes through these three points
p = polyfit(t_pts, y_pts, 2);    % p(1)*t^2 + p(2)*t + p(3)

% evaluate over full recording
F_model = polyval(p, t);

% align model mean to baseline
raw_base_mean   = mean(F_base);
model_base_mean = mean(F_model(baseline));
F_model_shifted = F_model - model_base_mean + raw_base_mean;

% ----- subtract bleaching model from signal -----
corrected_signal = F_sig - (F_model_shifted - raw_base_mean);

% ----- ΔF/F using corrected signal and baseline -----
F0 = mean(corrected_signal(baseline));
if F0 <= 0 || abs(F0) < 1e-6
    base_offset       = raw_base_mean;  % simple positive offset
    corrected_shifted = corrected_signal + base_offset;
    F0                = mean(corrected_shifted(baseline));
    dFF               = (corrected_shifted - F0) / F0;
else
    dFF = (corrected_signal - F0) / F0;
end

% ----- QC metrics -----
base_corr = corrected_signal(baseline);

% linear trend in baseline
p_lin = polyfit(t(baseline), base_corr, 1);
baseline_slope = p_lin(1);
base_corr_detrend = base_corr - polyval(p_lin, t(baseline));
baseline_std = std(base_corr_detrend);

% long-term drift (first vs last 600 s)
t_end_all = t(end);
early = t >= 0            & t <= 600;
late  = t >= (t_end_all-600)  & t <= t_end_all;
drift = mean(corrected_signal(late)) - mean(corrected_signal(early));

metrics.mode            = 'A';
metrics.F0              = F0;
metrics.raw_base_mean   = raw_base_mean;
metrics.model_base_mean = model_base_mean;
metrics.baseline_std    = baseline_std;
metrics.baseline_slope  = baseline_slope;
metrics.drift           = drift;
metrics.poly_coeffs     = p;   % store quadratic coefficients

% ----- optional diagnostic plot -----
figure;
subplot(4,1,1);
plot(t_base, F_base, 'k'); hold on;
plot(t_base, F_base_smooth, 'b');
plot(t, F_model_shifted, 'r','LineWidth',1.5);
legend('Baseline raw','Baseline smooth','Quadratic model');
title('Baseline region & quadratic model');
xlabel('Time (s)'); ylabel('F');

subplot(4,1,2);
plot(t, F_sig, 'k'); hold on;
plot(t, F_model_shifted, 'r');
legend('Signal','Bleach model');
title('Bleaching model over full recording');
xlabel('Time (s)'); ylabel('F');

subplot(4,1,3);
plot(t, corrected_signal, 'm');
xlabel('Time (s)'); ylabel('Corrected F');
title('Corrected signal (signal - bleach model)');

subplot(4,1,4);
plot(t, dFF, 'm');
xlabel('Time (s)'); ylabel('\DeltaF/F');
title('\DeltaF/F');
end