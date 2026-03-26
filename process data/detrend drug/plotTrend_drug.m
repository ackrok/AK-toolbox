function [ax, fig] = plotTrend_drug(out)
% Plot photometry signal trends and de-trended dF/F traces to evaluate
% varied methods for computing trend.
%
% Syntax:
%   [ax, fig] = plotTrend_drug(out);
% 
% Inputs:
%   out - output structure after running 'detrend_drug' fxn
% 
%       out.signal - photometry signal to detrend (vector)
%       out.Fs     - sampling rate, in Hz (scalar)
%       out.win    - time point, in seconds, that defines end of baseline (scalar)
%       out.trend  - table with varied baseline trends
%
% Outputs:
%   ax - handle to the axes created/selected by subplot
%           can use to further modify figure (hold on, set limits, add legends)
%
% Anya Krok, March 2026
%

%%
win = out.win;
signal = out.signal;
Fs = out.Fs;
trend = out.trend; 

%%
mask = 1 : win*Fs; % baseline mask
time = makeTime(length(signal), Fs); % time vector for plotting
ds = 1:Fs:numel(time); % downsampling to decrease sample rate for plotting

fig = figure; theme(fig, 'light');

%% (1) slow bleaching dynamics
ax(1) = subplot(2,2,1); hold on
plot(time(ds), signal(ds), 'k', 'DisplayName', 'raw signal');
plot(time(ds), trend.stretch(ds), 'g', 'LineWidth', 1.5, 'DisplayName', 'stretch');
plot(time(ds), trend.stitch(ds),  'b', 'LineWidth', 1.5, 'DisplayName', 'stitch');
plot(time(ds), trend.basesm(ds),  'c', 'LineWidth', 3, 'DisplayName', 'smooth baseline');
shadedband([win win+60], ylim, 'inj');
xlabel('time (s)'); ylabel('photometry (raw)');
title('baseline bleaching dynamics'); legend;

%% (2) double exponential model
ax(2) = subplot(2,2,2); hold on
plot(time(ds), signal(ds), 'k', 'DisplayName', 'raw signal');
plot(time(ds), trend.exp2stitch(ds), 'm', 'LineWidth', 2, 'DisplayName', 'exp2 on stitch');
plot(time(ds), trend.exp2base(ds), 'c--', 'LineWidth', 1.5, 'DisplayName', 'exp2 on baseline');
shadedband([win win+60], ylim, 'inj');
xlabel('time (s)'); ylabel('photometry (raw)');
title('double exponential model fits'); legend; 

%% compare model fits without data
% fig = figure; theme(fig,'light'); hold on
% plot(time(ds), trend.exp2stitch(ds), 'm', 'LineWidth', 1.5, 'DisplayName', 'exp2 on stitch');
% plot(time(ds), trend.stitch(ds),  'b', 'LineWidth', 1.5, 'DisplayName', 'stitched');
% shadedband([win win+60], ylim);
% xlabel('time (s)'); ylabel('photometry (raw)');
% title('comparing model fits');

%% (3, 4) compute and plot dF/F and z-score
var = trend.Properties.VariableNames; % trend labels
mat = table2array(trend); 
clr = {'g','b','m','c'};
ax(3) = subplot(2,2,3); hold on
ax(4) = subplot(2,2,4); hold on
for j = 1:4
    trend = mat(:,j); 
    [dff, z] = getdff_drug(out.y, trend, Fs, win);
    subplot(2,2,3); 
    plot(time(ds), dff(ds), clr{j}, 'DisplayName', var{j});
    subplot(2,2,4); 
    plot(time(ds), z(ds), clr{j}, 'DisplayName', var{j});
end
subplot(2,2,3); 
    shadedband([win win+60], ylim, 'inj');
    xlabel('time (s)'); ylabel('photometry (dF/F)'); legend; grid on
    title('comparing dF/F');
subplot(2,2,4); 
    shadedband([win win+60], ylim, 'inj');
    xlabel('time (s)'); ylabel('photometry (z-score)'); legend; grid on
    title('comparing z-score');

linkaxes(ax,'x');

end