%% Example trace of red signal detrending
% 
% JT030-260304: open field, saline injection
%
% pre-load variables
%   'signal' - raw photometry trace in red channel (red DA)
%   'Fs'     - sampling frequency in Hz (50 Hz)
%   'baseWin'- baseline range (600 sec)
%
%%
% a = 3; % ID within combined structure
% b = 2; % green or red channel

signal = comb(a).FP{b}; Fs = comb(a).Fs; 
baseWin = 600; 
time = makeTime(length(signal), Fs); % time vector for plotting
ds  = Fs; % downsampling by Fs for plotting
mask = 1 : baseWin*Fs; % baseline mask

%% process with dFF_drug
out = dFF_drug(signal, Fs, baseWin);

%%
fig = figure; theme(fig, 'light');

% (1) slow bleaching dynamics
subplot(2,3,1); hold on
plot(time(1:ds:end), out.sig(1:ds:end), 'k', ...
    'DisplayName', 'raw signal');
plot(time(1:ds:end), out.stretch(1:ds:end), 'g', 'LineWidth', 1.5, ...
    'DisplayName', 'stretch baseline');
plot(time(1:ds:end), out.stitch(1:ds:end),  'b', 'LineWidth', 1.5, ...
    'DisplayName', 'stitched');
plot(time(mask), out.basesm,  'c', 'LineWidth', 3, ...
    'DisplayName', 'smooth baseline');
shadedband([600 660], ylim);
legend
xlabel('time (s)'); ylabel('photometry (raw)');

% (2) double exponential model
subplot(2,3,2); hold on
plot(time(1:ds:end), signal(1:ds:end), 'k', ...
    'DisplayName', 'raw signal');
plot(time(1:ds:end), out.exp2stitch(1:ds:end), 'm', 'LineWidth', 2, ...
    'DisplayName', 'exp2 on stitch');
plot(time(1:ds:end), out.exp2base(1:ds:end), 'c--', 'LineWidth', 1.5, ...
    'DisplayName', 'exp2 on baseline');
shadedband([600 660], ylim);
legend
xlabel('time (s)'); ylabel('photometry (raw)');

% (3) compare models
subplot(2,3,3); hold on
plot(time(1:ds:end), out.exp2stitch(1:ds:end), 'm', 'LineWidth', 1.5, ...
    'DisplayName', 'exp2 on stitch');
plot(time(1:ds:end), out.stitch(1:ds:end),  'b', 'LineWidth', 1.5, ...
    'DisplayName', 'stitched');
shadedband([600 660], ylim);
legend
xlabel('time (s)'); ylabel('photometry (raw)');

% (4) dff
subplot(2,3,4); hold on
subplot(2,3,5); hold on
% (a) with stretch (prior model)
detrend = out.y - out.stretch;
F0 = mean(detrend(mask));
if F0 <= 0 || abs(F0) < 1e-3
    offset = mean(out.y(mask));   % use raw baseline level as DC offset
    detrend = detrend + offset;
    F0 = mean(detrend(mask)); 
end
sig_dff = (detrend - F0) ./ F0;
base_std = std(detrend(mask));
sig_z = (detrend - F0) ./ base_std;
subplot(2,3,4);
plot(time(1:ds:end), sig_dff(1:ds:end), 'g', 'DisplayName', 'stretch');
subplot(2,3,5);
plot(time(1:ds:end), sig_z(1:ds:end), 'g', 'DisplayName', 'stretch');

% (b) with stitch (new model)
detrend = out.y - out.stitch;
F0 = mean(detrend(mask));
if F0 <= 0 || abs(F0) < 1e-3
    offset = mean(out.y(mask));   % use raw baseline level as DC offset
    detrend = detrend + offset;
    F0 = mean(detrend(mask)); 
end
sig_dff = (detrend - F0) ./ F0;
base_std = std(detrend(mask));
sig_z = (detrend - F0) ./ base_std;
subplot(2,3,4);
plot(time(1:ds:end), sig_dff(1:ds:end), 'b', 'DisplayName', 'stretch');
subplot(2,3,5);
plot(time(1:ds:end), sig_z(1:ds:end), 'b', 'DisplayName', 'stretch');

% (c) with exp2stitch (newest model)
detrend = out.y - out.exp2stitch;
F0 = mean(detrend(mask));
if F0 <= 0 || abs(F0) < 1e-3
    offset = mean(out.y(mask));   % use raw baseline level as DC offset
    detrend = detrend + offset;
    F0 = mean(detrend(mask)); 
end
sig_dff = (detrend - F0) ./ F0;
base_std = std(detrend(mask));
sig_z = (detrend - F0) ./ base_std;
subplot(2,3,4);
plot(time(1:ds:end), sig_dff(1:ds:end), 'm--', 'DisplayName', 'stretch');
subplot(2,3,5);
plot(time(1:ds:end), sig_z(1:ds:end), 'm--', 'DisplayName', 'stretch');

subplot(2,3,4);
shadedband([600 660], ylim); grid on
legend; xlabel('time (s)'); ylabel('photometry (dF/F)');

subplot(2,3,5);
shadedband([600 660], ylim); grid on
legend; xlabel('time (s)'); ylabel('photometry (z-score)');

% title
subplot(2,3,1);
title(sprintf('%s-%s: %s',comb(a).mouse,comb(a).date,comb(a).FPnames{b}));