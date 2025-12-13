choice = menu('Cannula structure loaded into workspace?',...
    'YES','No, but will select from pre-saved cannula file','No, will load from raw data files');
switch choice
    case 2
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat'],'Select cannula file');
        load(fullfile(fPath,fName));
    case 3
        getCannula;
end
if ~exist('cannula','var')
    error('No cannula variable present in workspace');
end

%% INPUTS
nAn = length(cannula(1).s);
lickWithin = 0.25;      % CHANGE, lick within this window
winRew = [-1 2];        % CHANGE, window for aligning signal to rewarded trials, in seconds
winBase = [-4 -1];
winPkDA = [100 500];    % CHANGE, window for DA peak, in milliseconds
winTrACh = [300 700];   % CHANGE, window for ACh trough, in milliseconds
xlbl = 'reward';
winAcc = [-1 1];        % CHANGE, window for aligning signal to acceleration, in seconds
NumStd = 2;             % CHANGE, for immobility trough/peak analysis
clr = {'k','r','b','c','m'};
choice = menu('Choose infusion to compare to saline', {cannula.inf}); % select infusions to compare from those present in cannula structure
choice = [1, choice]; % compare to saline (should be first in list)

%% Align photometry to reward
for z = 1:length(cannula)
    beh = cannula(z).s;
    r = cannula(z).win;
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        if isempty(beh(x).reward) 
            cannula(z).a_rew{x,1} = nan(151,1);
            cannula(z).a_rew{x,2} = nan(151,1);
            continue; end
        rew = beh(x).reward./Fs; % reward delivery times, in seconds
        lick = beh(x).lick./Fs; % lick times, in seconds
        [rewYes, ~, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
        %rewYes: indices of rewarded trials (animal licked within specified window)
        ev = rew(rewYes); % delivery times for rewarded trials, in seconds
        ev = ev(ev > r(x,1)*60 & ev < r(x,2)*60); % extract rewarded trials during post-infusion window
        
%         % % aligning photometry to FIRST lick
%         bin = 1/1000;
%         peth = getClusterPETH(lickNew, ev, bin, [0 1]); % PETH: lick aligned to rewarded trials in 1 ms bins
%         cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
%         [~, lickFirst] = max(cts~=0, [], 1); % Find first non-zero index for each trial
%         lickFirst = lickFirst(:).*bin + ev(:); % First lick after delivery for rewarded trials
%         ev = lickFirst; 
%         xlbl = '1st lick'; winPkDA = [0 400]; winTrACh = [200 600];
        
        for y = 1:length(beh(x).FP)
            sig = beh(x).FP{y}; % signal that will be aligned to event times
            sig = sig - nanmean(sig); % subtract mean of trace to center on zero
            [mat, time] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]);
            [matBase] = getSTA(sig, ev, Fs, [winBase(1), winBase(end)]);
            matBase = nanmean(matBase,1); % average across baseline window
            cannula(z).a_rew{x,y} = mat - matBase; % save into structure
            
            sig = sig(randperm(length(sig),length(sig))); % random permutation of photometry signal
            [matShuff] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]);
            [matBase] = getSTA(sig, ev, Fs, [winBase(1), winBase(end)]);
            cannula(z).a_shuff{x,y} = matShuff - matBase; % save into structure
        end
    end
end

%% PLOT average reward responses
if ~isfield(cannula(1),'a_rew')
    error('No field a_rew - must align photometry to reward delivery before proceeding.');
end
fig = figure; fig.Position(3) = 1000;
for z = 1:2
    a_rew = cannula(choice(z)).a_rew; % extract reward-aligned data from structure
    for y = 1:size(a_rew,2)
        a_mat = []; % initialize temporary output matrix
        for x = 1:size(a_rew,1)
            a_mat(:,x) = nanmean(a_rew{x,y},2); % extract reward-aligned data from structure and average across all rewards
        end
        subplot(1,2,y); hold on
        plot(time, a_mat, 'Color', [char2rgb(clr{choice(z)}), 0.2]);
        % shadederrbar(time, nanmean(a_mat,2), SEM(a_mat,2), clr{choice(z)}); % plot average across trials
        ylabel(sprintf('%s amplitude (dF/F)',cannula(1).s(1).FPnames{y}));
        xlabel(sprintf('time to %s (s)',xlbl));
        title(sprintf('%s (n = %d): %s vs %s',...
            cannula(1).s(1).FPnames{y},...
            size(a_mat(:,~isnan(a_mat(1,:))),2),...
            cannula(choice(1)).inf,cannula(choice(2)).inf)); 
        axis square
    end
end
movegui(gcf,'center');

%% PLOT example reward responses
x = 2; % example
if ~isfield(cannula(1),'a_rew')
    error('No field a_rew - must align photometry to reward delivery before proceeding.');
end
fig = figure; fig.Position(3) = 1000;
for z = 1:2
    a_rew = cannula(choice(z)).a_rew; % extract reward-aligned data from structure
    for y = 1:size(a_rew,2)
        subplot(1,2,y);
        if z == 1 % plot shuffled saline photometry aligned to reward
            shuff = cannula(z).a_shuff{x,y};
            shadederrbar(time, nanmean(nanmean(shuff,2)).*ones(length(time),1), nanmean(SEM(shuff,2)).*ones(length(time),1), 'b');
        end
        shadederrbar(time, nanmean(a_rew{x,y},2), SEM(a_rew{x,y},2), clr{choice(z)}); % plot average across trials
        ylabel(sprintf('%s amplitude (dF/F)',cannula(1).s(1).FPnames{y}));
        xlabel(sprintf('time to %s (s)',xlbl));
        title(sprintf('%s (%s): %s vs %s',...
            cannula(1).s(1).FPnames{y},...
            strtok(cannula(1).s(x).rec,'-'),...
            cannula(choice(1)).inf,cannula(choice(2)).inf)); 
        axis square
    end
end
movegui(gcf,'center');

%% Find AMPLITUDE and LATENCY of reward response
if ~isfield(cannula(1),'a_rew')
    error('No field a_rew - must align photometry to reward delivery before proceeding.');
end
lag = {}; val = {}; % initialize output
for z = 1:2
    a_rew = cannula(choice(z)).a_rew; % extract reward-aligned data from structure
    Fs = cannula(choice(z)).s(1).Fs; % sampling frequency
    y = 1; % analyzing ACh reward response
    win = winTrACh; % range for ACh trough, in milliseconds
    r = (win/(1000/Fs) + find(time == 0)); 
    for x = 1:size(a_rew,1)
        pull = a_rew{x,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
        pull = pull(r(1):r(2),:); % extract only segment of time specified as window for reward-related trough
        [a,b] = min(pull); % find local MINIMUM within range
        c = time(b + r(1) - 1); % convert index to seconds
        c(isnan(a)) = nan; 
        val{y}(x,z) = nanmean(a); % save amplitude of minima
        lag{y}(x,z) = nanmean(c).*1000; % save lag at minima, in milliseconds
    end
    val{y}(val{y} == 0) = nan; lag{y}(lag{y} == 0) = nan;
    
    y = 2; % analyzing DA reward response
    win = winPkDA; % range for ACh trough, in milliseconds
    r = (win/(1000/Fs) + find(time == 0)); 
    for x = 1:size(a_rew,1)
        pull = a_rew{x,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
        pull = pull(r(1):r(2),:); % extract only segment of time specified as window for reward-related trough
        [a,b] = max(pull); % find local MAXIMUM within range
        c = time(b + r(1) - 1); % convert index to seconds
        c(isnan(a)) = nan; 
        val{y}(x,z) = nanmean(a); % save amplitude of minima
        lag{y}(x,z) = nanmean(c).*1000; % save lag at minima, in milliseconds
    end
    val{y}(val{y} == 0) = nan; lag{y}(lag{y} == 0) = nan;
end

%% PLOT AMPLITUDE + LATENCY OF REWARD RESPONSE
fig = figure; fig.Position(3) = 1000;
for y = 1:2
    subplot(1,2,y); hold on
    for z = 1:2
        for x = 1:nAn
            plot(lag{y}(x,:), val{y}(x,:), '-', 'Color', [0 0 0 0.1],'HandleVisibility','off');
        end
        plot(lag{y}(:,z), val{y}(:,z), '.', 'MarkerSize', 20, 'Color', clr{choice(z)}); % plot individual data points for each mouse
        errorbar(nanmean(lag{y}(:,z)), nanmean(val{y}(:,z)), ... % plot average across all mice
        SEM(val{y}(:,z),1), SEM(val{y}(:,z),1), ... % error bars vertical for amplitude
        SEM(lag{y}(:,z),1), SEM(lag{y}(:,z),1), ... % error bars horizontal for latency
        '.', 'MarkerSize', 20, 'Color', clr{choice(z)},'HandleVisibility','off');
        xlabel(sprintf('time to %s (s)',xlbl)); xlim([win(1) win(2)]);
        legend({cannula.inf},'Location','southwest');
    end
end
y = 1; subplot(1,2,y); % ACh reward response subplot
    xlabel('time to reward (s)'); xlim([winTrACh(1) winTrACh(2)]);
    ylabel('ACh trough amplitude (%dF/F)'); ylim([-4 0]); yticks([-8:1:0]);
    [~,p] = ttest(lag{y}(:,1),lag{y}(:,2)); % statistical test: paired t-test
    [~,p(2)] = ttest(val{y}(:,1),val{y}(:,2)); % statistical test: paired t-test
    title(sprintf('%s vs %s (lag %1.2f, val %1.2f)',cannula(choice(1)).inf,cannula(choice(2)).inf,p(1),p(2)));
    axis square
y = 2; subplot(1,2,y); % DA reward response subplot
    xlabel('time to reward (s)'); xlim([winPkDA(1) winPkDA(2)]);
    ylabel('DA peak amplitude (%dF/F)'); ylim([0 15]);
    [~,p] = ttest(lag{y}(:,1),lag{y}(:,2)); % statistical test: paired t-test
    [~,p(2)] = ttest(val{y}(:,1),val{y}(:,2)); % statistical test: paired t-test
    title(sprintf('%s vs %s (lag %1.2f, val %1.2f)',cannula(choice(1)).inf,cannula(choice(2)).inf,p(1),p(2)));
    axis square
movegui(gcf,'center');

% subplot(1,3,3); hold on
% y = 1;
% plot(val{y}','.-r', 'MarkerSize', 20);
% errorbar([0.75 2.25],nanmean(val{y}),SEM(val{y},1),'.k', 'MarkerSize', 20);
% xlim([0.5 2.5]); xticks([1 2]); xticklabels({cannula.inf});
% ylim([-4 0]); yticks([-4:0]); ylabel('ACh trough amplitude (%dF/F)');
% [~,p] = ttest(val{y}(:,1),val{y}(:,2)); % statistical test: paired t-test
% title(sprintf('%s vs %s (val p = %1.2f)',cannula(choice(1)).inf,cannula(choice(2)).inf,p));
% axis square