function [out] = detrendFP_drug(signal, Fs, win)
% Detrend photometry signal using pre-drug baseline.
%
%   [out] = detrendFP_drug(signal, Fs, win)
%
% Description: This code will detrend photometry signal using values from 
% specified baseline period (eg, pre drug injection) with a variety of
% methods.
%
% Input:
%    signal - photometry signal to detrend (vector)
%    Fs     - sampling rate, in Hz (scalar)
%    win    - time point, in seconds, that defines end of baseline (scalar)
%       --> baseline will be signal(1:winBase*Fs);
%
% Output:
%    out - structure with variables including
%       out.signal - input signal (vector)
%       out.stretch, stitch, exp2stitch, exp2base - varied trends fits (vectors)
%       out.lbl - name of trend used to compute dFF, z-score (char array)
%       out.dff - processed photometry signal, units deltaF/F (vector)
%       out.z   - processed photometry signal, units z-score (vector)
%
% Author: Anya Krok, March 2026
% Last Updated: March 24 2026

%%
plotY  = false; % plotting logical, default false
signal = signal(:); % ensure column vector for faster computation
nSig   = numel(signal); % number of samples
time   = makeTime(nSig, Fs); % create time vector

%% define parameters
sm = 100; % sliding window length for smoothing baseline, in seconds

%% initial processing of signal
signal(1:round(Fs/2)) = signal(round(Fs/2)); % remove artifact of 
% resampling in first half second of session

% low-pass filter data
cutoff = 15; % cut-off frequency for filter
order  = 8; % order of the filter
y_full = filterFP(signal,Fs,cutoff,order,'lowpass');
y_full = y_full(:); 

%% define baseline
mask   = 1 : win*Fs; % create baseline mask
t_base = time(mask(:));
y_base = y_full(mask(:)); 
nBase  = numel(y_base);

%% create artificial bleaching curve using baseline
% (1) smoothing baseline to keep only slow bleaching dynamics
basesm = movmean(y_base, round(sm*Fs)); 
% (2) stretch (aka resample) smooth baseline to full length
stretch = resample(basesm, nSig, nBase); % stretched baseline trend
stretch = stretch(:);
stretch = movmean(stretch, sm*Fs); % smooth edges

%% stitched baseline trend
stitch = zeros(nSig,1); % initialize
stitch(mask) = basesm; % use smooth non-stretched baseline for baseline period
% At baseline end, stretch has offset from baseline due to resample
% process. In order to stitch without jump will need to subtract offset.
maskBaseEnd = nBase-Fs : nBase; % mask for last second of baseline
offset = mean(stretch(maskBaseEnd)) - mean(basesm(maskBaseEnd));
stretch_offset = stretch - offset;
stitch(nBase+1:end) = stretch_offset(nBase+1:end); % add to output
% Smooth vector to remove discontinuities, residual steps
bw = 10; % transition window, in seconds
j0 = nBase;
j1 = min(nSig, nBase + bw);
tau = (0:(j1-j0)) / max(1,(j1-j0));
wblend = 0.5*(1 + cos(pi * tau)); % 1 -> 0
for k = 0:(j1-j0); idx = j0 + k; a = wblend(k+1);
    stitch(idx) = a * basesm(min(idx, nBase)) + (1 - a) * stretch_offset(idx);
end
stitch = movmean(stitch, sm*Fs); % smoothing again

%% fit double exponential model to stitch
y_toFit = stitch;
fo = fitoptions('exp2',...
    'Lower', [-Inf, -Inf, -Inf, -Inf], ... % constrain so exponents are non-positive (force decays)
    'Upper', [Inf, 0, Inf, 0]);
f = fit(time, y_toFit, 'exp2', fo); % can compute coefficients over full session because y-input is full stitched trend
exp2stitch = feval(f, time); % evaluate over full session

%% fit double exponential model to baseline values
y_toFit = basesm;
fo = fitoptions('exp2', ...
    'Lower', [-Inf, -Inf, -Inf, -Inf], ... % constrain so exponents are non-positive (force decays)
    'Upper', [Inf, 0, Inf, 0], ...
    'StartPoint', [max(y_toFit)-median(y_toFit), -1/win, 0.5*(max(y_toFit)-median(y_toFit)), -0.1/win]);
f = fit(t_base, y_toFit, 'exp2', fo); % compute coefficients only over baseline period
exp2base = feval(f, time); % evaluate over full session

%% compute dF/F and z-score
% select fit to detrend signal
trend = exp2stitch; lbl = 'exp2stitch'; % CHANGE
% trend = stretch; lbl = 'stretch'; % CHANGE
% trend = stitch; lbl = 'stitch'; % CHANGE
% trend = exp2base; lbl = 'exp2base'; % CHANGE
%
detrend = y_full(:) - trend(:);
%
% compute dF/F
F0 = mean(detrend(mask));
if F0 <= 1e-3 % if F0 is too small or negative will add offset
    offset = mean(y_base); % use average of baseline as DC offset
    detrend = detrend + offset;
    F0 = F0 + offset;
end
sig_dff = (detrend - F0) ./ F0;
%
% compute z-score
sigma = std(detrend(mask));
sig_z = (detrend - F0) ./ sigma;

%% output 
out = struct;
out.signal = signal;  % input signal
out.Fs  = Fs;      % input sampling frequency
out.win = win; % baseline window
out.y   = y_full;  % signal after low-pass filtering
out.dff = sig_dff; % units: dF/F
out.z   = sig_z;   % units: z-score
out.lbl = lbl;     % label of trend utilized for dff, z-score
out.basesm = basesm;   % smooth baseline
out.stretch = stretch; % stretch trend
out.stitch = stitch;   % stitched trend (smooth baseline with offset stretch)
out.exp2stitch = exp2stitch; % double exp model fit to stitch trend
out.exp2base = exp2base;     % double exp model fit to baseline values only

%% PLOT
if plotY
    fig = figure; theme(fig, 'light');
    sp(1) = subplot(3,1,1);
    plot(time, signal, 'g'); hold on;
    plot(time, trend, 'r', 'LineWidth',2);
    xline(win,'LineWidth',2);
    xlabel('Time (s)'); ylabel('F');
    legend('rawFP', 'baseline');
    title('raw photometry and baseline fit');
    
    sp(2) = subplot(3,1,2); hold on
    plot(time, sig_dff, 'm');
    xline(win,'LineWidth',2);
    xlabel('Time (s)'); ylabel('\DeltaF/F');
    title('\DeltaF/F');
    
    sp(3) = subplot(3,1,3); hold on
    plot(time, sig_z, 'm');
    xline(win,'LineWidth',2);
    xlabel('Time (s)'); ylabel('z-score');
    title('z-score');
    
    linkaxes(sp,'x');
end

%% Plot to inspect
% fig = figure; theme(fig, 'light'); hold on;
% plot(time, signal, 'k', 'DisplayName','signal');
% plot(t_base, basesm, 'c--', 'DisplayName','baseline smoothed');
% plot(time, stitch, 'b', 'DisplayName','stitch');
% plot(time, exp2stitch, 'r', 'DisplayName','exp2stitch');
% plot(time, exp2base, 'g', 'DisplayName','exp2base');
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
