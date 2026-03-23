function [detrend, finalTrend, fitted_exp2] = analyzeSignalBiexp(signal, Fs, baseWin)
% analyzeSignalBiexp  Fit data-driven biexponential baseline and detrend signal
% 
% Usage:
%   [detrend, fittedTrend] = analyzeSignalBiexp(signal, Fs, baseWin)
%
% Inputs:
%   signal   - column vector (n x 1)
%   time     - column vector (n x 1) in seconds
%   Fs       - sampling rate (Hz)
%   baseWin  - scalar, baseline window length in seconds
%
% Optional inputs:
%   'Plot'        - true/false (default false)
%
% Outputs:
%   detrend      - detrended signal (n x 1)
%   fittedTrend  - fitted baseline trend (n x 1)
%   fitted_exp2  - fitted baseline trend with exp2

%%
doPlot = 'false';
signal = signal(:);
time = makeTime(numel(signal),Fs);

% baseline indices
maskN = max(3, round(baseWin*Fs));
mask = 1:maskN;
t_base = time(mask) - time(mask(1));
y_base = signal(mask);

% light smoothing to remove high-frequency noise (helps local slope estimation)
y_sm = movmean(y_base, max(1,round(0.2*Fs)));

% sliding local log-fit to estimate local decay rates k(t)
win_dur = min(baseWin/4, 2);                 % seconds, tweakable
winN = max(5, round(win_dur*Fs));
stepN = max(1, round(0.2*Fs));
starts = 1:stepN:(numel(y_sm)-winN+1);
kvals = nan(numel(starts),1);

for ii = 1:numel(starts)
    idx = starts(ii):(starts(ii)+winN-1);
    yw = y_sm(idx);
    tw = time(idx) - time(idx(1));
    % estimate local offset as median of window tail
    Cw = median(yw(max(1,end-round(0.1*winN)+1):end));
    yC = yw - Cw;
    if all(yC > 0) && numel(tw) >= 3
        pfit = polyfit(tw, log(yC), 1);
        kvals(ii) = -pfit(1);
    end
end

kvals = kvals(kvals>0 & isfinite(kvals));
% fallback if insufficient estimates
if numel(kvals) < 6
    k_mean = max(1e-6, nanmedian(kvals));
    k_fast = k_mean;
    k_slow = k_mean;
else
    % cluster into two groups
    try
        [~, C] = kmeans(kvals(:), 2, 'Replicates',5, 'MaxIter',200);
        % cluster centroids: larger k -> faster decay
        k_fast = max(C);    % faster is larger k
        k_slow = min(C);
        % ensure k_fast >= k_slow
        if k_fast < k_slow
            [k_fast, k_slow] = deal(k_slow, k_fast);
        end
    catch
        % if kmeans fails
        k_mean = median(kvals);
        k_fast = k_mean; k_slow = k_mean;
    end
end

tau_fast = 1/max(k_fast,1e-12);
tau_slow = 1/max(k_slow,1e-12);

% build biexponential fittype and sensible start points
ft = fittype('A*exp(-k1*x) + B*exp(-k2*x) + C', ...
             'independent','x','coefficients',{'A','k1','B','k2','C'});
opts = fitoptions(ft);
opts.Lower = [-Inf, 0, -Inf, 0, -Inf];
% amplitudes: split baseline drop into A0,B0; k guesses from taus
A0 = (y_base(1)-median(y_base(end-max(2,round(0.1*maskN))+1:end))) * 0.6;
B0 = (y_base(1)-median(y_base(end-max(2,round(0.1*maskN))+1:end))) * 0.4;
k1_0 = 1 / (max(0.1, 0.5*tau_fast)); % heuristic
k2_0 = 1 / (max(0.1, 2*tau_slow));
C0 = median(y_base(end-max(5,round(0.05*maskN))+1:end));

opts.StartPoint = [A0, k1_0, B0, k2_0, C0];
opts.Robust = 'LAR';

% fit, with fallback to monoexponential if biexp fails
fittedTrend_baseline = nan(maskN,1);
fitObj = []; gof = [];
try
    [fobj,gof,~] = fit(t_base, y_base, ft, opts);
    fittedTrend_baseline = fobj(t_base);
    fitObj = fobj; gof = gof;
catch
    % try monoexponential A*exp(-k*x)+C
    ft1 = fittype('A*exp(-k*x) + C','independent','x','coefficients',{'A','k','C'});
    opts1 = fitoptions(ft1);
    opts1.Lower = [-Inf, 0, -Inf];
    opts1.StartPoint = [A0+B0, 1/(max(0.1, mean([tau_fast,tau_slow]))), C0];
    opts1.Robust = 'LAR';
    try
        [fobj,gof,~] = fit(t_base, y_base, ft1, opts1);
        fittedTrend_baseline = fobj(t_base);
        fitObj = fobj; gof = gof;
    catch
        % final fallback: smooth baseline and use stretched version
        fittedTrend_baseline = movmean(y_base, max(1,round(1*Fs)));
        fitObj = [];
        gof = [];
    end
end

% evaluate trend across full recording (time relative to baseline start)
t_full = time - time(mask(1));
if ~isempty(fitObj)
    fittedTrend = fitObj(t_full);
else
    % stretch smoothed baseline to full length
    fittedTrend = interp1(t_base + time(mask(1)), fittedTrend_baseline, time, 'pchip', 'extrap');
end
% optional mild smoothing to remove tiny wiggles
fittedTrend = movmean(fittedTrend, max(1,round(2*Fs)));
finalTrend = fittedTrend;

% % robust linear alignment of fittedTrend to full signal
% % uses only times where fittedTrend and signal are finite
% ok = isfinite(fittedTrend) & isfinite(signal);
% y = signal(ok);
% x = fittedTrend(ok);
% 
% % use robust fit to reduce outlier influence
% try
%     % fit y = a*x + b robustly
%     w = robustfit(x, y);   % robustfit returns [b; a]
%     b = w(1); a = w(2);
% catch
%     % fallback to ordinary least squares
%     p = polyfit(x, y, 1);  % p(1) = a, p(2) = b
%     a = p(1); b = p(2);
% end
% 
% % constrain a to be near 1 if you expect only offset differences:
% a = max(0.8, min(1.2, a)); % optional bound
% fittedTrend_aligned = a * fittedTrend + b;
% finalTrend = fittedTrend_aligned;

%%
% produce detrended output
detrend = signal - finalTrend;

% outputs
fitInfo.fitObj = fitObj;
fitInfo.gof = gof;
fitInfo.tau_fast = tau_fast;
fitInfo.tau_slow = tau_slow;
fitInfo.k_fast = k_fast;
fitInfo.k_slow = k_slow;
fitInfo.startGuess = opts.StartPoint;

% comparison with built-in 'exp2' model
% Use 'exp2' model; set bounds to encourage decay (B,D <= 0)
% 'exp2' uses A*exp(B*x)+C*exp(D*x)
opts2 = fitoptions('Method','NonlinearLeastSquares');
f_exp2 = fit(t_base, y_base, 'exp2', opts2);
% set generous bounds but prefer negative exponents (decay)
opts2.Lower = [-Inf, -Inf, -Inf, -Inf];
opts2.Upper = [Inf, 0, Inf, 0];   % require B<=0, D<=0 (decaying)
opts2.StartPoint = []; % let fit choose
opts2.Robust = 'LAR';
fitted_exp2 = f_exp2(t_full);
fitted_exp2 = movmean(fitted_exp2, max(1, round(2*Fs)));

% optional plot for diagnostics
    % fig = figure; theme(fig,'light');
    % subplot(2,1,1);
    % plot(time, signal, 'k'); hold on;
    % plot(time, fittedTrend, 'r','LineWidth',1.5);
    % plot(time, fittedTrend_aligned, 'c','LineWidth',1.5);
    % plot(time, fitted_exp2, 'm','LineWidth',1.5);
    % xlim([time(1), time(end)]);
    % xlabel('Time (s)'); ylabel('Signal');
    % legend('raw','fitted trend','exp2');
    % title('Signal and fitted baseline');
    % subplot(2,1,2);
    % plot(t_base + time(mask(1)), y_base, '.b'); hold on;
    % plot(t_base + time(mask(1)), fittedTrend_baseline, '-r','LineWidth',1.2);
    % xlabel('Time (s)'); ylabel('Baseline segment');
    % legend('baseline data','baseline fit');

end