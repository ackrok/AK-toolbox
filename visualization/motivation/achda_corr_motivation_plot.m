%% PLOT N = X mice
fig = figure; fig.Position(3) = 1375; 
clr = {'r','k';'g','k';'b','k'};
lbl = {'immobility','locomotion','reward'};
for y = 1:3
    subplot(1,3,y); hold on
    for z = 1:2
        a = corr_an{y,z}; a = a - nanmean(a(find(lags./Fs == -2):find(lags./Fs == -0.5),:));
        shadederrbar(lags/Fs, nanmean(a,2), SEM(corr_an{y,z},2), clr{y,z}); % EARLY 
    end
    plot([0 0],[-1 0.5],'--k');
    a = nanmean(corr_shuff{2,y},2); b = nanmean(corr_shuff{3,y},2);
    a = a - nanmean(a(find(lags./Fs == -2):find(lags./Fs == -1),:));
    shadederrbar(lags/Fs, a, b, 'k'); hold on % SHUFFLE
    xlabel('Lag (s)'); xlim([-0.5 0.5]);
    ylabel('Correlation Coefficient'); ylim([-1 0.5]); yticks([-1:0.25:1]);
    title(sprintf('%s (n = %d mice)',lbl{y}, nAn));
    axis square; set(gca,'TickDir','out');
    legend({'early','late','shuff'});
end
movegui(gcf,'center');

%% PLOT STATS - early vs late for each state
fig = figure; fig.Position([3 4]) = [1000 620];
clr = {'r','k';'g','k';'b','k'};
lbl = {'imm','loc','rew'}; lbl2 = {'early','late'};
a = 0.8; b = 1.2; r1 = a + (b-a).*rand(nAn,1); % jitter for plotting early data points
a = 1.8; b = 2.2; r2 = a + (b-a).*rand(nAn,1); % jitter for plotting late data points

for y = 1:3
    subplot(2,3,y); hold on
    pull = [min_val(:,y), min_val_late(:,y)]; % plotting 1, plotting 2
    plot([r1';r2'], pull', '.k', 'MarkerSize', 20);
    errorbar([1,2], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{y,1});
    xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl2); 
    ylabel('Max Coefficient'); ylim([-1 0]); yticks([-1:0.25:0]);
    [~,p] = ttest(pull(:,1),pull(:,2));
    title(sprintf('%s (p = %1.3f)',lbl{y},p)); axis square
    
    subplot(2,3,y+3); hold on
    pull = [min_lag(:,y), min_lag_late(:,y)].*1000; % plotting 1, plotting 2
    plot([r1';r2'], pull', '.k', 'MarkerSize', 20);
    errorbar([1,2], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{y,1});
    xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl2); 
    ylabel('Latency (ms)'); ylim([-300 0]); yticks([-500:100:0]);
    [~,p] = ttest(pull(:,1),pull(:,2));
    title(sprintf('%s (p = %1.3f)',lbl{y},p)); axis square
end
movegui(gcf,'center');

%% PLOT STATS - early by state, late by state
fig = figure; fig.Position([3 4]) = [650 620];
clr = {'r','k';'g','k';'b','k'};
lbl = {'imm','loc','rew'}; lbl2 = {'early','late'};
a = 0.8; b = 1.2; r1 = a + (b-a).*rand(nAn,1); % jitter for plotting immobility
a = 1.8; b = 2.2; r2 = a + (b-a).*rand(nAn,1); % jitter for plotting locomotion
a = 2.8; b = 3.2; r3 = a + (b-a).*rand(nAn,1); % jitter for plotting reward

for y = 1:2
    subplot(2,2,y); hold on
    switch y; case 1; pull = min_val; case 2; pull = min_val_late; end
    plot([r1';r2';r3'], pull', '.k', 'MarkerSize', 20);
    errorbar([1 2 3], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{y,1});
    xlim([0.5 3.5]); xticks([1 2 3]); xticklabels(lbl); 
    ylabel('Max Coefficient'); ylim([-1 0]); yticks([-1:0.25:0]);
    [p, ~, stats] = anova1(pull,[],'off'); % anova
    c = multcompare(stats,'display','off');
    title(sprintf('%s (anova = %1.3f)',lbl2{y},p)); axis square
    
    subplot(2,2,y+2); hold on
    switch y; case 1; pull = min_lag; case 2; pull = min_lag_late; end
    pull = pull.*1000;
    plot([r1';r2';r3'], pull', '.k', 'MarkerSize', 20);
    errorbar([1 2 3], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{y,1});
    xlim([0.5 3.5]); xticks([1 2 3]); xticklabels(lbl); 
    ylabel('Latency (ms)'); ylim([-300 0]); yticks([-500:100:0]);
    [p, ~, stats] = anova1(pull,[],'off'); % anova
    c = multcompare(stats,'display','off');
    title(sprintf('%s (anova = %1.3f)',lbl2{y},p)); axis square
end
movegui(gcf,'center');

%% PLOT STATS -- within animal comparison of correlation coefficients and latency to maximum
fig = figure; fig.Position([3 4]) = [1000 620];
clr = {'r','k';'g','k';'b','k'};
lbl = {'imm','loc','rew'}; lbl2 = {'early','late'};

for y = 1:3
    subplot(2,3,y); hold on
    pull = [min_val(:,y), min_val_late(:,y)]; % plotting 1, plotting 2
    a = nanmean(corr_an{y,3},2);
    a = a - nanmean(a(find(lags/Fs == -5):find(lags/Fs == -1),:));
    shadederrbar([0.5;2.5], zeros(2,1), nanmean(nanmean(corr_an{y,4}-corr_an{y,3},2)).*ones(2,1), 'k'); % shuffle
    plot(pull', '.-k', 'MarkerSize', 20); % maximum coefficient for each animal, connecting early-late
    errorbar([0.8 2.2], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{y,1}); % mean +/- SEM
    xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl2); 
    ylabel('Max Coefficient'); ylim([-1 0.25]); yticks([-1:0.25:0]);
    [~,p] = ttest(pull(:,1),pull(:,2));
    title(sprintf('%s (p = %1.3f)',lbl{y},p)); axis square
    
    subplot(2,3,y+3); hold on
    pull = [min_lag(:,y), min_lag_late(:,y)].*1000; % plotting 1, plotting 2
    plot(pull', '.-k', 'MarkerSize', 20);
    errorbar([0.8 2.2], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{y,1});
    xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl2); 
    ylabel('Latency (ms)'); ylim([-300 0]); yticks([-500:100:0]);
    [~,p] = ttest(pull(:,1),pull(:,2));
    title(sprintf('%s (p = %1.3f)',lbl{y},p)); axis square
end
movegui(gcf,'center');
