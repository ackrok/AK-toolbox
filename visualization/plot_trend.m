% Plot photometry signal trends and de-trended dF/F traces to evaluate
% varied methods for computing trend.
% Will run detrendFP_drug.
% 
% Inputs:
%    signal - photometry signal to detrend (vector)
%    Fs     - sampling rate, in Hz (scalar)
%    win    - time point, in seconds, that defines end of baseline (scalar)
%
% Anya Krok, March 2026
%
%%
% a = 3; % mouse ID
% b = 2; % green or red channel
% signal = comb(a).FP{b}; Fs = comb(a).Fs; 
win = 600;

if exist('out','var') && isstruct(out) && isfield(out,'trend')
    % proceed with script
else
    assert(exist('signal','var') & exist('Fs','var'), ...
        'Error: Load signal into workspace.');
    out = detrendFP_drug(signal, Fs, win); % RUN
    fprintf('Done: detrendFP_drug using signal, injTime %ds).\n',win);
end

%%
mask = 1 : out.win*Fs; % baseline mask
time = makeTime(length(out.signal), out.Fs); % time vector for plotting
ds = 1:Fs:numel(time); % downsampling to decrease sample rate for plotting

fig = figure; theme(fig, 'light');

% (1) slow bleaching dynamics
ax(1) = subplot(2,2,1); hold on
plot(time(ds), out.signal(ds), 'k', 'DisplayName', 'raw signal');
plot(time(ds), out.trend.stretch(ds), 'g', 'LineWidth', 1.5, 'DisplayName', 'stretch');
plot(time(ds), out.trend.stitch(ds),  'b', 'LineWidth', 1.5, 'DisplayName', 'stitch');
plot(time(ds), out.trend.basesm(ds),  'c', 'LineWidth', 3, 'DisplayName', 'smooth baseline');
shadedband([win win+60], ylim);
xlabel('time (s)'); ylabel('photometry (raw)');
title('baseline bleaching dynamics'); legend;

% (2) double exponential model
ax(2) = subplot(2,2,2); hold on
plot(time(ds), out.signal(ds), 'k', 'DisplayName', 'raw signal');
plot(time(ds), out.trend.exp2stitch(ds), 'm', 'LineWidth', 2, 'DisplayName', 'exp2 on stitch');
plot(time(ds), out.trend.exp2base(ds), 'c--', 'LineWidth', 1.5, 'DisplayName', 'exp2 on baseline');
shadedband([win win+60], ylim);
xlabel('time (s)'); ylabel('photometry (raw)');
title('double exponential model fits'); legend; 

% (3) compare model fits without data
ax(3) = subplot(2,2,3); hold on
plot(time(ds), out.trend.exp2stitch(ds), 'm', 'LineWidth', 1.5, 'DisplayName', 'exp2 on stitch');
plot(time(ds), out.trend.stitch(ds),  'b', 'LineWidth', 1.5, 'DisplayName', 'stitched');
shadedband([win win+60], ylim);
xlabel('time (s)'); ylabel('photometry (raw)');
title('comparing model fits');

% (4) compute and plot dF/F
var = out.trend.Properties.VariableNames; % trend labels
mat = table2array(out.trend); 
clr = {'g','b','m','c'};
ax(4) = subplot(2,2,4); hold on
for j = 1:4
    detrend = out.y - mat(:,j); % subtract trend
    F0 = mean(detrend(mask));
    if F0 <= 1e-3
        offset  = mean(out.y(mask)); % DC offset
        detrend = detrend + offset;
        F0 = F0 + offset;
    end
    sig_dff = (detrend - F0) ./ F0;
    % sigma = std(detrend(mask));
    % sig_z = (detrend - F0) ./ sigma;
    subplot(2,2,4); 
    plot(time(ds), sig_dff(ds), clr{j}, 'DisplayName', var{j});
end
shadedband([win win+60], ylim);
xlabel('time (s)'); ylabel('photometry (dF/F)'); legend; grid on
title('comparing dF/F');
linkaxes(ax,'x');