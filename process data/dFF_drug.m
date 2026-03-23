function [out] = dFF_drug(signal, Fs, baseWin)
% Process photometry signal using predrug baseline.
%
%   [out] = dFF_drug(signal, Fs, baseWin)
%
% Description: This code will baseline adjust the photometry signal using
% baseline period, prior to drug administration
%
% Input:
%    signal - photometry signal to baseline (column vector)
%    Fs - sampling rate, in Hz (scalar)
%    winBase - baseline range from 1:X in seconds (scalar)
%
% Output:
%    out - structure with output variables including
%       out.dff - processed photometry signal, as deltaF/F
%       out.z   - processed photometry signal, z-scored
%    fittedTrend - fitted baseline trend
%
% Author: Anya Krok, March 2026
%%
plotY  = false; % plotting logical, default false
signal = signal(:); % ensure column vector for faster computation
nSig   = numel(signal); % number of samples
time   = makeTime(nSig, Fs); % create time vector

% remove artifact of resampling in first half second of session
signal(1:round(Fs/2)) = signal(round(Fs/2));

% low-pass filter data
cutoff = 15; % cut-off frequency for filter
order = 8; % order of the filter
y_full = filterFP(signal,Fs,cutoff,order,'lowpass');
y_full = y_full(:); 

% define baseline
mask = 1 : baseWin*Fs; % create baseline mask
t_base = time(mask(:));
y_base = signal(mask(:)); 
nBase = numel(y_base);

%% create artificial bleaching curve using baseline
sm = 100; % smooth baseline to keep only slow bleaching dynamics
basesm = movmean(y_base, round(sm*Fs)); 
% stretch (resample) smooth baseline to full length
stretch = resample(basesm, nSig, nBase); % stretched baseline trend
stretch = stretch(:);
stretch = movmean(stretch, sm*Fs); % smooth edges

%% combo stitch baseline and stretch at transition point
% March 21 2026
maskBaseEnd = nBase-Fs : nBase; % mask for section of baseline
offset = mean(stretch(maskBaseEnd)) - mean(basesm(maskBaseEnd)); % offset for stretch above baseline end
adj = stretch - offset;
%
stitch = zeros(nSig,1); % initialize
stitch(mask) = basesm; % use smooth non-stretched baseline for baseline period
stitch(nBase+1:end) = adj(nBase+1:end); % remainder of signal
%
bw = 10; % transition window to smooth residual step, in seconds
j0 = nBase;
j1 = min(nSig, nBase + bw);
tau = (0:(j1-j0)) / max(1,(j1-j0));
wblend = 0.5*(1 + cos(pi * tau)); % 1 -> 0
for k = 0:(j1-j0); idx = j0 + k; a = wblend(k+1);
    stitch(idx) = a * basesm(min(idx, nBase)) + (1 - a) * adj(idx);
end
stitch = movmean(stitch, sm*Fs); % 1 s smoothing to remove tiny discontinuities

%% fit a double exponential model to the stitch
% March 23 2026
y_forexp2 = stitch;
fo = fitoptions('exp2',...
    'Lower', [-Inf, -Inf, -Inf, -Inf], ... % constrain so exponents are non-positive (force decays)
    'Upper', [Inf, 0, Inf, 0]);
f = fit(time, y_forexp2, 'exp2', fo); % perform fit on stitch
exp2stitch = feval(f, time);

%% double exponential model
% March 20 2026
y_forexp2 = movmean(y_base, round(50*Fs)); 
fo = fitoptions('exp2', ...
    'Lower', [-Inf, -Inf, -Inf, -Inf], ... % constrain so exponents are non-positive (force decays)
    'Upper', [Inf, 0, Inf, 0], ...
    'StartPoint', [max(y_base)-median(y_base), -1/baseWin, 0.5*(max(y_base)-median(y_base)), -0.1/baseWin]);
f = fit(t_base, y_forexp2, 'exp2', fo);
exp2base = f(time);

%% compute dF/F
trend = exp2stitch; lbl = 'exp2stitch'; % CHANGE
% trend = stretch; lbl = 'stretch'; % CHANGE
% trend = stitch; lbl = 'stitch'; % CHANGE
% trend = exp2base; lbl = 'exp2base'; % CHANGE

detrend = y_full(:) - trend(:); 
F0 = mean(detrend(mask));
% If F0 is too close to zero or negative, add a constant offset
if F0 <= 0 || abs(F0) < 1e-3
    offset = mean(y_full(mask));   % use raw baseline level as DC offset
    detrend = detrend + offset;
    F0 = mean(detrend(mask)); 
end
sig_dff = (detrend - F0) ./ F0;

%% z-score
base_std = std(detrend(mask));
sig_z = (detrend - F0) ./ base_std;

%% output 
out = struct;
out.sig = signal;  % input signal
out.Fs  = Fs;      % input sampling frequency
out.win = baseWin; % baseline window
out.y   = y_full;  % signal after low-pass filtering
out.dff = sig_dff; % units: dF/F
out.z   = sig_z;   % units: z-score
out.trend = lbl;   % label of trend utilized for dff, z-score
out.basesm = basesm; % smooth baseline
out.stretch = stretch; % stretch trend
out.stitch = stitch;   % stitch trend of smooth baseline and offset stretch
out.exp2stitch = exp2stitch; % double exp model fit from stitch trend
out.exp2base = exp2base;    % double exp model fit from baseline

%% PLOT
if plotY
    fig = figure; theme(fig, 'light');
    sp(1) = subplot(3,1,1);
    plot(time, signal, 'g'); hold on;
    plot(time, trend, 'r', 'LineWidth',2);
    xline(baseWin,'LineWidth',2);
    xlabel('Time (s)'); ylabel('F');
    legend('rawFP', 'baseline');
    title('raw photometry and baseline fit');
    
    sp(2) = subplot(3,1,2); hold on
    plot(time, sig_dff, 'm');
    xline(baseWin,'LineWidth',2);
    xlabel('Time (s)'); ylabel('\DeltaF/F');
    title('\DeltaF/F');
    
    sp(3) = subplot(3,1,3); hold on
    plot(time, sig_z, 'm');
    xline(baseWin,'LineWidth',2);
    xlabel('Time (s)'); ylabel('z-score');
    title('z-score');
    
    linkaxes(sp,'x');
end

%% Plot to inspect
% fig = figure; theme(fig, 'light');
% plot(time, signal, 'k', 'DisplayName','signal'); hold on;
% plot(time(1:nBase), base, 'bo-', 'DisplayName','baseline smoothed (no stretch)');
% plot(time, exp2, 'g--', 'LineWidth',2,'DisplayName','exp2'); 
% plot(time, stretch, 'r--', 'DisplayName','stretched baseline (full)');
% plot(time, trend, 'm-', 'LineWidth',1.6, 'DisplayName','combined fittedTrend');
% xlim([min(time) max(time)]);
% legend('Location','best'); xlabel('Time'); ylabel('Signal / trend');
% title('Baseline: smooth (baseline) stitched with stretched baseline (post)');

%% OLDER ATTEMPTS

%% complex exponential model
% March 20 2026
% base = exp2(mask);
% stretch = resample(basesm, nSig, nBase); % stretched trend
% stretch = stretch(:);
% stretch = movmean(stretch, sm*Fs); % smooth edges
% maskBaseEnd = nBase-Fs : nBase; % mask for section of baseline
% offset = mean(stretch(maskBaseEnd)) - mean(base(maskBaseEnd)); % offset for stretch above baseline end
% adj = stretch - offset;
% %
% fittedTrend = zeros(nSig,1); % initialize
% fittedTrend(mask) = basesm; % use smooth non-stretched baseline for baseline period
% fittedTrend(nBase+1:end) = adj(nBase+1:end); % remainder of signal
% %
% bw = 10; % transition window to smooth residual step, in seconds
% j0 = nBase;
% j1 = min(nSig, nBase + bw);
% tau = (0:(j1-j0)) / max(1,(j1-j0));
% wblend = 0.5*(1 + cos(pi * tau)); % 1 -> 0
% for k = 0:(j1-j0); idx = j0 + k; a = wblend(k+1);
%     fittedTrend(idx) = a * basesm(min(idx, nBase)) + (1 - a) * adj(idx);
% end
% fittedTrend = movmean(fittedTrend, sm*Fs); % 1 s smoothing to remove tiny discontinuities
% detrend = y_full(:) - fittedTrend(:);

%% double exponential model
% March 20 2026
% [detrend, trend, fitted_exp2] = analyzeSignalBiexp(signal, Fs, baseWin);

%% monoexponential model
% % define monoexponential model: A*exp(-k*t) + C
% ft = fittype('A*exp(-k*x) + C','independent','x','coefficients',{'A','k','C'});
% opts = fitoptions(ft);
% opts.Lower = [-Inf, 0, -Inf];        % enforce nonnegative decay rate k>=0
% opts.StartPoint = [y_base(1)-y_base(end), 1/baseWin, mean(y_base(end-1:end))];
% % fit (robust option helps if outliers)
% [fh, ~] = fit(t_base, y_base, ft, opts);
% % evaluate fit across full time and smooth short-scale residuals
% sm = 5; % in seconds
% trend = fh(time(:));            % fitted trend over full recording
% trend = movmean(trend, max(1,round(sm*Fs))); % Xs smoothing
% % subtract to detrend
% detrend = y_full(:) - trend;        % baseline-subtracted

%%
% %Baseline signal
% % use baseline FP to identify baseline of specified window using a
% % moving window and finding values within a specified percentile
% interpType = 'linear'; % 'linear' 'spline' 
% fitType = 'interp';  % Fit method 'interp' , 'exp' , 'line'
% basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline
% %Note: Lower prc are used because the mean of signal is not true baseline
% winSize = 10; % Window size for baselining in seconds
% winOv = 0; %Window overlap size in seconds
% [~, baseline] = baselineFP (signal,interpType,fitType,basePrc,winSize,winOv,Fs);
% %
% signal = signal - baseline; % subtract baesline
% 
% %Baseline of predrug portion of signal
% predrug = signal(mask); % extract portion of signal specified by window
% basePrc = 50; % Percentile value from 1 - 100 to use when finding baseline
% [~, predrug_baseline] = baselineFP (predrug,interpType,fitType,basePrc,winSize,winOv,Fs);
% predrug_mu = mean(predrug_baseline);
%
% % Compute dF/F of entire signal
% dFF = (signal - predrug_mu) ./ predrug_mu;
% dFF = dFF * 100; 

end
