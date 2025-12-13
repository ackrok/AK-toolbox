%% LOAD DATA
% Comparing immobility ACh pause peak statistics
% for 2 different sets of recordings

beh1 = extractBeh; % Load PREVIOUS reward recordings
beh2 = extractBeh; % Load NEW reward recordings
compare = struct;
compare(1).s = beh1; compare(1).lbl = 'pre';
compare(2).s = beh2; compare(2).lbl = 'post';

%% INPUTS
nAn = length(compare(1).s);
NumStd = 2; % CHANGE, for immobility trough/peak analysis

%% Immobility pause peak
for y = 1:length(compare)
    beh = compare(y).s; % Extract behavior
    switch y
        case 1
            [amp, dur, freq, thres] = getImmPausePeak(beh, NumStd); % Immobility pause peak statistics
        case 2
            % [amp, dur, freq, ~] = getImmPausePeak(beh, thres);
            [amp, dur, freq, thres] = getImmPausePeak(beh, NumStd);
    end
    Fs = beh(1).Fs;
    freq(isnan(freq)) = nan;
    dur(isnan(dur)) = nan; dur = dur.*(1000/Fs); % Adjust from samples to ms
    amp(isnan(amp)) = nan; amp = abs(amp); % Adjust to be absolute amplitude
    rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
    uni = unique(rec); nAn = length(uni); % number of unique animals
    plotme = cell(3,1);
    for x = 1:nAn
        ii = find(strcmp(rec,uni{x}));
        plotme{1}(x,:) = nanmean(freq(ii,:),1); % average across multiple recordings from same mouse
        plotme{2}(x,:) = nanmean(dur(ii,:),1); % average across multiple recordings from same mouse
        plotme{3}(x,:) = nanmean(amp(ii,:),1); % average across multiple recordings from same mouse
    end
    compare(y).out = plotme;
    compare(y).lbl1 = {'trough','peak'};
    compare(y).lbl2 = {'frequency','duration','amplitude'};
end
fprintf('Done. \n'); 

%%
fig = figure; fig.Position([3 4]) = [1375 800];
lbl1 = compare(1).lbl1;
lbl2 = compare(1).lbl2;

nAn = length(beh);
lims1 = [0 2; 0 500; 0 8];
ticks1 = {[0:0.5:2], [0:100:500], [0:2:8]};
for jj = 1:2 % Iterate over pause and peaks
    for ii = 1:length(plotme)
        subplot(2,3,ii + (jj-1)*3); hold on
        a = [compare(1).out{ii}(:,jj), compare(2).out{ii}(:,jj)];
        plot(a', '.-', 'Color', [0 0 0 0.1], 'MarkerSize', 20);
        errorbar([0.75 2.25], nanmean(a), SEM(a,1), '.r', 'MarkerSize', 20);
        xlim([0.5 2.5]); xticks([1 2]); xticklabels({compare.lbl});
        ylabel(sprintf('%s',lbl2{ii})); ylim(lims1(ii,:)); yticks(ticks1{ii});
        [~,p] = ttest(a(:,1),a(:,2));
        title(sprintf('%s-%s (p = %1.3f)',lbl1{jj},lbl2{ii},p));
        axis square; set(gca,'TickDir','out');
    end
end
movegui(gcf,'center');
