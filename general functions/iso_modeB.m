function [corrected_signal, dFF, metrics] = iso_modeB(signal, fs, modelType)
% iso_modeB
%   Mode B: Use first 10 min of the *signal itself* to fit a bleaching model.
%   - Clamp first 2 s
%   - Take 0–600 s baseline from signal
%   - Lightly smooth baseline (to reduce noise)
%   - Fit a mono-exponential: a*exp(-b*t) + c on baseline
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
    error('iso_modeB: no samples in baseline window 0–%d s.', baseline_sec);
end

t_base = t(baseline);
F_base = F_sig(baseline);

% ----- lightly smooth baseline to reduce noise -----
smooth_win_sec     = 10;                        % 10 s smoothing (tune as needed)
smooth_win_samples = max(1, round(smooth_win_sec * fs));
F_base_smooth      = movmean(F_base, smooth_win_samples);

% ----- fit bleaching model on baseline -----
% normalize time on baseline to [0,1] to improve numerical stability
tb0     = t_base(1);
tb1     = t_base(end);
tb_norm = (t_base - tb0) / (tb1 - tb0);   % 0..1 over baseline

switch lower(modelType)
    case 'exp1'  % mono-exponential: a*exp(-b*x) + c
        ft  = fittype('a*exp(-b*x) + c', 'independent','x', 'dependent','y');
        opts = fitoptions(ft);

        a0 = F_base_smooth(1) - F_base_smooth(end);
        b0 = 1;
        c0 = F_base_smooth(end);

        opts.StartPoint = [a0, b0, c0];
        opts.Lower      = [-Inf, 0, -Inf];
        opts.Upper      = [ Inf, 10, Inf];

    case 'exp2'  % double-exponential: a1*exp(-b1*x) + a2*exp(-b2*x) + c
        ft  = fittype('a1*exp(-b1*x) + a2*exp(-b2*x) + c', ...
                      'independent','x', 'dependent','y');
        opts = fitoptions(ft);

        % crude initial guesses
        a1_0 = (F_base_smooth(1) - F_base_smooth(end))/2;
        a2_0 = a1_0;
        b1_0 = 1;
        b2_0 = 0.1;
        c0   = F_base_smooth(end);

        opts.StartPoint = [a1_0, b1_0, a2_0, b2_0, c0];
        opts.Lower      = [-Inf, 0, -Inf, 0, -Inf];
        opts.Upper      = [ Inf, 10,  Inf, 10,  Inf];

    case 'poly1' % linear: p1*x + p2
        ft  = fittype('p1*x + p2', 'independent','x', 'dependent','y');
        opts = fitoptions(ft);

        p_lin = polyfit(tb_norm, F_base_smooth, 1);
        opts.StartPoint = p_lin;

    case 'poly2' % quadratic: p1*x^2 + p2*x + p3
        ft  = fittype('p1*x^2 + p2*x + p3', 'independent','x', 'dependent','y');
        opts = fitoptions(ft);

        p_quad = polyfit(tb_norm, F_base_smooth, 2);
        opts.StartPoint = p_quad;

     case 'poly3' % cubic: p1*x^3 + p2*x^2 + p3*x + p4
        ft  = fittype('p1*x^3 + p2*x^2 + p3*x + p4', ...
                      'independent','x', 'dependent','y');
        opts = fitoptions(ft);

        % Use polyfit to get a reasonable initial guess
        p_cubic = polyfit(tb_norm, F_base_smooth, 3);  % returns [p1 p2 p3 p4]
        opts.StartPoint = p_cubic;

    otherwise
        error('Unknown modelType "%s". Use exp1, exp2, poly1, poly2, ...', modelType);
end

exp_fit = fit(tb_norm, F_base_smooth, ft, opts);

% ----- evaluate bleaching model over full recording -----
t0     = t(1);
t1     = t(end);
t_norm = (t - t0) / (t1 - t0);       % 0..1 over full recording

F_model = exp_fit(t_norm);           % same length as F_sig

% optional: align model mean to signal mean on baseline
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
t_end = t(end);
early = t >= 0            & t <= 600;
late  = t >= (t_end-600)  & t <= t_end;
drift = mean(corrected_signal(late)) - mean(corrected_signal(early));

metrics.mode           = 'B';
metrics.F0             = F0;
metrics.raw_base_mean  = raw_base_mean;
metrics.model_base_mean = model_base_mean;
metrics.baseline_std   = baseline_std;
metrics.baseline_slope = baseline_slope;
metrics.drift          = drift;
metrics.exp_fit        = exp_fit;     % in case you want to inspect parameters

% ----- optional diagnostic plot -----
figure;
subplot(4,1,1);
plot(t_base, F_base, 'k'); hold on;
plot(t_base, F_base_smooth, 'b');
plot(t_base, exp_fit(tb_norm), 'r','LineWidth',1.5);
legend('Baseline raw','Baseline smooth','Exp fit (baseline)');
title('Baseline region');
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
