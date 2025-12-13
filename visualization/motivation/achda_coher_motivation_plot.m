%% EARLY vs LATE coherence - plot average across animals

fig = figure; fig.Position([3 4]) = [1000 620];
clr = {'r','k';'g','k';'b','k'}; lbl = {'immobility','locomotion','reward'};
for y = 1:3 % iterate across behavioral states
    subplot(2,3,y); hold on
    for z = 1:2 % iterate over early / late
        shadederrbar(f, nanmean(coher_an{y,z},2), SEM(coher_an{y,z},2), clr{y,z}); end
    plot([0.5 0.5],[0 1],'--k'); plot([4 4],[0 1],'--k');
    shadederrbar(f, nanmean(coher_shuff{2},2), nanmean(coher_shuff{3},2) - nanmean(coher_shuff{2},2), 'k');
    % legend({'early','','late','','0.5Hz','4Hz'});
    ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
    xlabel('frequency');
    title(sprintf('%s',lbl{y})); axis square

    subplot(2,3,y+3); hold on
    for z = 1:2 % iterate over early / late
        shadederrbar(f, nanmean(rad2deg(-phase_an{y,z}),2), SEM(rad2deg(-phase_an{y,z}),2), clr{y,z}); end
    plot([0.5 0.5],[-180 180],'--k'); plot([4 4],[-180 180],'--k');
    shadederrbar(f, nanmean(rad2deg(phase_shuff{2}),2), nanmean(rad2deg(phase_shuff{3}),2) - nanmean(rad2deg(phase_shuff{2}),2), 'k');
    % legend({'early','','late','','0.5Hz','4Hz'});
    ylabel('degrees'); ylim([-180 180]); yticks([-180:90:180]);
    title(sprintf('%s phase',lbl{y})); axis square
end
movegui(gcf,'center');

%% EARLY vs LATE coherence for each STATE - statistics
fig = figure; fig.Position([3 4]) = [1000 620];
clr = {'r','k';'g','k';'b','k'};
lbl = {'imm','loc','rew'}; lbl2 = {'early','late'};
r = [6:42]; % [~,r2] = min(abs(f - 2)); % range for 0.5-4Hz
a = 0.8; b = 1.2; r1 = a + (b-a).*rand(nAn,1); % jitter for plotting early data points
a = 1.8; b = 2.2; r2 = a + (b-a).*rand(nAn,1); % jitter for plotting late data points
for y = 1:3 % iterate across behavioral states
    coher_avg = []; phase_avg = [];
    for z = 1:2 % iterate over early / late
        % coher_avg(:,z) = coher_an{y,z}(r2,:); phase_avg(:,z) = phase_an{y,z}(r2,:); % value at 2Hz
        % coher_avg(:,z) = nanmean(coher_an{y,z}(r,:)); phase_avg(:,z) = nanmean(phase_an{y,z}(r,:)); % average within frequency band
        coher_avg(:,z) = median(coher_an{y,z}(r,:)); 
        phase_avg(:,z) = median(phase_an{y,z}(r,:)); % median within frequency band
    end
    subplot(2,3,y); hold on 
    pull = coher_avg;
    plot([r1';r2'], pull', '.k', 'MarkerSize', 20);
    errorbar([1,2], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{y,1});
    xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl2); 
    ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
    [~,p] = ttest(pull(:,1),pull(:,2));
    title(sprintf('%s coher (p = %1.3f)',lbl{y},p)); axis square

    subplot(2,3,y+3); hold on
    pull = rad2deg(phase_avg);
    plot([r1';r2'], pull', '.k', 'MarkerSize', 20);
    errorbar([1 2], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{y,1});
    xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl2); 
    ylabel('degrees'); ylim([0 180]); yticks([-180:90:180]);
    [~,p] = ttest(pull(:,1),pull(:,2));
    title(sprintf('%s phase (p = %1.3f)',lbl{y},p)); axis square
end
movegui(gcf,'center');

%% EARLY coher by STATE, LATE coher by STATE - statistics
fig = figure; fig.Position([3 4]) = [650 620];
clr = {'r','k';'g','k';'b','k'};
lbl = {'imm','loc','rew'}; lbl2 = {'early','late'};
a = 0.8; b = 1.2; r1 = a + (b-a).*rand(nAn,1); % jitter for plotting immobility
a = 1.8; b = 2.2; r2 = a + (b-a).*rand(nAn,1); % jitter for plotting locomotion
a = 2.8; b = 3.2; r3 = a + (b-a).*rand(nAn,1); % jitter for plotting reward

for z = 1:2
    coher_avg = []; phase_avg = [];
    for y = 1:3 % iterate over early / late
        % coher_avg(:,z) = coher_an{y,z}(r2,:); phase_avg(:,z) = phase_an{y,z}(r2,:); % value at 2Hz
        % coher_avg(:,z) = nanmean(coher_an{y,z}(r,:)); phase_avg(:,z) = nanmean(phase_an{y,z}(r,:)); % average within frequency band
        coher_avg(:,y) = median(coher_an{y,z}(r,:)); 
        phase_avg(:,y) = median(phase_an{y,z}(r,:)); % median within frequency band
    end
    
    subplot(2,2,z); hold on
    pull = coher_avg;
    plot([r1';r2';r3'], pull', '.k', 'MarkerSize', 20);
    errorbar([1 2 3], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{z,1});
    xlim([0.5 3.5]); xticks([1 2 3]); xticklabels(lbl); 
    ylabel('coherence'); ylim([0 1]); yticks([0:0.2:1]);
    [p, ~, stats] = anova1(pull,[],'off'); % anova
    c = multcompare(stats,'display','off');
    title(sprintf('%s (anova = %1.3f)',lbl2{z},p)); axis square
    
    subplot(2,2,z+2); hold on
    pull = rad2deg(phase_avg);
    plot([r1';r2';r3'], pull', '.k', 'MarkerSize', 20);
    errorbar([1 2 3], nanmean(pull), SEM(pull,1), '.', 'MarkerSize', 20, 'Color', clr{z,1});
    xlim([0.5 3.5]); xticks([1 2 3]); xticklabels(lbl); 
    ylabel('degrees'); ylim([0 180]); yticks([-180:90:180]);
    [p, ~, stats] = anova1(pull,[],'off'); % anova
    c = multcompare(stats,'display','off');
    title(sprintf('%s (anova = %1.3f)',lbl2{z},p)); axis square
end
movegui(gcf,'center');