%% alternative
% high-pass filter to remove slow decay but preserve drug dynamics
% period 10-30 mins for drug response that evolves over seconds-a few min
% period 30-60 mins if drug effect is slower, evolves over several minutes
% if drug onset amplitude is reduced or distorted, try LONGER period
% for ketamine data, period of ~50 mins works well

b = 1;
for a = 1:3

    Fs = comb(a).Fs;  % sampling rate, per seconds
    signal = comb(a).nbFP{b};
    time   = makeTime(numel(signal), Fs); % create time vector
    
    % params
    filtOrder = 2;    % filter order
    period = 50;      % period in minutes (change)
    win = 10;         % window in minutes (change)
    
    % high-pass filter
    fc = 1/period;    % cycles per minute
    cutoff = fc/((Fs*60)/2); % normalized cutoff for butter (Nyquist in samples/min)
    [coeffB,coeffA] = butter(filtOrder, cutoff, 'high'); % high-pass filter
    y_hp = filtfilt(coeffB,coeffA,signal);
    
    % fit double exponential model to baseline values
    win = win * 60; % baseline window, in seconds
    mask   = 1 : win*Fs; mask = mask(:); % baseline mask
    y_base = y_hp(mask);
    fopts = fitoptions('exp2', 'Lower', [-Inf, -Inf, -Inf, -Inf], 'Upper', [Inf, 0, Inf, 0]);
    %    'StartPoint', [max(y_toFit)-median(y_toFit), -1/win, 0.5*(max(y_toFit)-median(y_toFit)), -0.1/win]);
    f = fit(time(mask), y_base, 'exp2', fopts); % compute coefficients over baseline period
    exp2base = feval(f, time); % evaluate over full session
    
    % detrend and compute dF/F and z-score
    trend = exp2base;
    detrend = y_hp(:) - trend(:);
    detrend = detrend + abs(min(detrend)); % offset so no negative values
    F0 = mean(detrend(mask)); % compute dF/F:
    if F0 <= 1e-3 % if F0 is too small or negative will add offset
        offset  = mean(y_base); % compute a DC offset
        detrend = detrend + offset;
        F0 = F0 + offset;
    end
    dff = (detrend - F0) ./ F0;
    sigma = std(detrend(mask));
    sig_z = (detrend - F0) ./ sigma;
    
    fig = figure; theme(fig,'light'); 
    subplot(4,1,1); hold on; grid on
    title(sprintf('%s-%s demodulated signal',comb(a).mouse,comb(a).date));
    plot(time, signal,'k','LineWidth',2);
    subplot(4,1,2); hold on; grid on
    title(sprintf('high-pass filter signal (period: %d min)',period));
    plot(time, y_hp,'b','DisplayName',sprintf('%d',period));
    plot(time, exp2base,'c','LineWidth',2);
    subplot(4,1,3); hold on; grid on
    title(sprintf('dF/F (baseline: %d min)',win/60));
    plot(time, dff,'m');
    subplot(4,1,4); hold on; grid on
    title('z-score');
    plot(time, sig_z,'g');
    
    comb(a).dff = dff;
    comb(a).z = sig_z;
end

%% 
% % testing different periods:
% fig = figure; theme(fig,'light'); 
% subplot(2,1,1);
% plot(signal,'k','LineWidth',2);
% subplot(2,1,2); hold on
% 
% opts = [5 10 20 30 60]; % period options
% for a = 1:5
%     period = opts(a); % period, in minutes
%     fc = 1/period;    % cycles per minute
%     Fs = 50;          % sampling rate, per seconds
%     cutoff = fc/((Fs*60)/2); % normalized cutoff for butter (Nyquist in samples/min)
%     [coeffB,coeffA] = butter(filtOrder, cutoff, 'high'); % high-pass filter
%     y_hp = filtfilt(coeffB,coeffA,signal);
%     plot(movmean(y_hp,50),'DisplayName',sprintf('%d',period));
% end
% legend
