function [out] = detrendFP_drug(signal, Fs, win, varargin)
% This function will detrend photometry signal using a trend curve fitted
% to the data in specified baseline window (eg, pre drug injection).
%
% Syntax:
%   [out] = detrendFP_drug(signal, Fs, win);
%   [out] = detrendFP_drug(signal, Fs, win, trendLbl);
%
% Inputs:
%    signal - photometry signal to detrend (vector)
%    Fs     - sampling rate, in Hz (scalar)
%    win    - time point, in seconds, that defines end of baseline (scalar)
%       --> baseline will be signal(1:winBase*Fs);
%    trendLbl - (optional) label for trend to use to compute dF/F, z-score
%       character array or string array
%       options: 'stretch', 'stitch', 'exp2stitch', 'exp2base'
%
% Outputs:
%    out - structure with variables including
%       out.signal - input signal (vector)
%       out.lbl - label of trend used to compute dFF, z-score (char array)
%       out.dff - processed photometry signal, units deltaF/F (vector)
%       out.z   - processed photometry signal, units z-score  (vector)
%       out.trend - table with fitted trends
%
% Filtering Options: default low-pass filter cutoff 15, order 8
% Smoothing Options: default 100 s
% Trend Options: use trendLbl if is input, or default 'exp2stitch'
%
% Author: Anya Krok, March 2026
% Last Updated: March 24 2026

%% parameters
cutoff = 15; % filter: cut-off frequency,   default 15 Hz
order  = 8;  % filter: order of the filter, default 8
sm = 100;    % smooth baseline: sliding window length, default 100 s
trendlbl = 'exp2stitch'; % trend to use, default 'exp2stitch'
if ~isempty(varargin), trendlbl = varargin{1}; end % override if provided
plotY  = false; % plotting logical, default false

%% 
signal = signal(:); % ensure column vector for faster computation
nSig   = numel(signal); % number of samples
time   = makeTime(nSig, Fs); % create time vector

%% filter
% cutoff = 15; % cut-off frequency for filter
% order  = 8; % order of the filter
y_full = filterFP(signal,Fs,cutoff,order,'lowpass');
y_full = y_full; 

%% define baseline
mask   = 1 : win*Fs; mask = mask(:); % baseline mask
tbase = time  (mask);
ybase = y_full(mask); 
nBase  = numel (ybase);

%% bleaching curve using smooth baseline
% (1) compute smooth baseline, to retain slow bleaching dynamics
% Create interpolated line using the bottom percentile of values
basePrc = 1; winSize = 1; winOv = 0; interpType = 'linear'; fitType = 'interp';
[~,basebot] = baselineFP(ybase,interpType,fitType,basePrc,winSize,winOv,Fs);
basesm = movmean(basebot, round(sm*Fs)); % then smooth

% (2) expand slow baseline trend to match length of signal with 'resample'
basesm(1:round(Fs/2)) = basesm(round(Fs/2)); % flatten first half-second
% Why? resample may create edge artifacts at start of sequence, replacing
% first 0.5s with constant value will blunt transient prior to resampling.
stretch = resample(basesm, nSig, nBase); 
% 'resample' function will expand sample count, but preserve low-frequency
% content and applies anti-aliasing filter
stretch = stretch(:);
stretch = movmean(stretch, sm*Fs); % smooth and ripples introduced by resample and smooth edges

%% stitched bleaching curve
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
fopts = fitoptions('exp2',...
    'Lower', [-Inf, -Inf, -Inf, -Inf], ... % constrain so exponents are non-positive
    'Upper', [Inf, 0, Inf, 0]);
f = fit(time, y_toFit, 'exp2', fopts); % can compute coefficients over full session (y-input is full stitched trend)
exp2stitch = feval(f, time); % evaluate over full session

%% fit double exponential model to baseline values
y_toFit = basesm;
fopts = fitoptions('exp2', 'Lower', [-Inf, -Inf, -Inf, -Inf], 'Upper', [Inf, 0, Inf, 0]);
%    'StartPoint', [max(y_toFit)-median(y_toFit), -1/win, 0.5*(max(y_toFit)-median(y_toFit)), -0.1/win]);
f = fit(tbase, y_toFit, 'exp2', fopts); % compute coefficients over baseline period
exp2base = feval(f, time); % evaluate over full session

%% compute dF/F and z-score
switch trendlbl
    case 'exp2stitch', trend = exp2stitch; % default unless other input is provided 
    case 'stretch', trend = stretch;
    case 'stitch', trend = stitch; 
    case 'exp2base', trend = exp2base;
end
detrend = y_full(:) - trend(:);
% compute dF/F:
F0 = mean(detrend(mask));
if F0 <= 1e-3 % if F0 is too small or negative will add offset
    offset  = mean(ybase); % compute a DC offset
    detrend = detrend + offset;
    F0 = F0 + offset;
end
sig_dff = (detrend - F0) ./ F0;
% compute z-score:
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
out.lbl = trendlbl;     % label of trend utilized for dff, z-score
% %
varNames = {'stretch','stitch','exp2stitch','exp2base','basesm'};
C   = {stretch, stitch, exp2stitch, exp2base, basesm}; % load into cell array
% %
C   = cellfun(@(x) x(:), C, 'UniformOutput', false);   % ensure column vectors
len = cellfun(@numel, C); maxlen = max(len); 
if any(len ~= maxlen) % pad shorter cells with NaN if needed
    C = cellfun(@(x) [x; nan(maxlen - numel(x), 1)], C, 'UniformOutput', false);
end
T = table(C{:}, 'VariableNames', varNames);
out.trend = T; % load table into output structure

% Table includes:
% out.trend.stretch; % stretch trend
% out.trend.stitchh; % stitched trend (smooth baseline with offset stretch)
% out.trend.exp2stitch; % double exp model fit to stitch trend
% out.trend.exp2base;   % double exp model fit to baseline values only
% out.trend.basesm;  % smooth baseline

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

end
