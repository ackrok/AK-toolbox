% check on signal
%%
win = 800;
nFPchan = 1;
opts = {'exp2stitch','exp2base','stitch','stretch'};

for b = 1:nFPchan
    for a = 1:length(comb)
        if isfield(comb,'nbFP'); signal = comb(a).nbFP{b};
            else, signal = comb(a).FP{b}; 
        end
        Fs = comb(a).Fs;
        out = detrend_drug(signal, Fs, win);
        [ax, fig] = plotTrend_drug(out);
        ID = [comb(a).mouse,'-',comb(a).date];
        choice = menu(sprintf('Select trend for dF/F: %s',ID),opts);
        chosen = opts{choice};
        dff = getdff_drug(out.y, out.trend.(chosen), Fs, win);
        comb(a).dff(:,b) = dff;
        comb(a).trend{b} = chosen; 
        close(fig);
    end
end

%% alternative
% high-pass filter to remove slow decay but preserve drug dynamics
% period 10-30 mins for drug response that evolves over seconds-a few min
% period 30-60 mins if drug effect is slower, evolves over several minutes
% if drug onset amplitude is reduced or distorted, try LONGER period
% for ketamine data, period of ~50 mins works well

% a = 1; b = 1;
% win = 12*60; % seconds for window
% 
% Fs = comb(a).Fs;
% signal = comb(peraiod).nbFP{b};
% 
% fig = figure; theme(fig,'light'); 
% subplot(2,1,1);
% plot(signal,'k','LineWidth',2);
% subplot(2,1,2); hold on
filtOrder = 2; % filter order
opts = [5 10 20 30 60]; % period options
for a = 1:5
period = opts(a); % period, in minutes
fc = 1/period;    % cycles per minute
Fs = 50;          % sampling rate, per seconds
cutoff = fc/((Fs*60)/2); % normalized cutoff for butter (Nyquist in samples/min)
[coeffB,coeffA] = butter(filtOrder, cutoff, 'high'); % high-pass filter
y_hp = filtfilt(coeffB,coeffA,signal);
plot(movmean(y_hp,50),'DisplayName',sprintf('%d',period));
end
legend

win = 800;
nFPchan = 1;
opts = {'exp2stitch','exp2base','stitch','stretch'};

for b = 1:nFPchan
    for a = 1:length(comb)
        if isfield(comb,'nbFP'); signal = comb(a).nbFP{b};
            else, signal = comb(a).FP{b}; 
        end
        Fs = comb(a).Fs;
        out = detrend_drug(signal, Fs, win);
        [ax, fig] = plotTrend_drug(out);
        ID = [comb(a).mouse,'-',comb(a).date];
        choice = menu(sprintf('Select trend for dF/F: %s',ID),opts);
        chosen = opts{choice};
        dff = getdff_drug(out.y, out.trend.(chosen), Fs, win);
        comb(a).dff(:,b) = dff;
        comb(a).trend{b} = chosen; 
        close(fig);
    end
end