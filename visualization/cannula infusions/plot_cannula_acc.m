loaded = menu('Already loaded cannula into workspace?','yes','no');
switch loaded
    case 2
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat']); % Select .mat files you want to load
        load(fullfile(fPath,fName));
end
if ~exist('cannula','var')
    getCannula
end

%% INPUTS
nAn = length(cannula(1).s);
lickWithin = 0.25;      % CHANGE, lick within this window
winRew = [-1 2];        % CHANGE, window for aligning signal to rewarded trials, in seconds
winBase = [-4 -1];
winAccPk = [-100 300; 100 500];    % CHANGE, window for ACh/DA peak, in milliseconds
xlbl = 'acceleration';
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
        acc = getAcc(beh(x).vel); % Extract acceleration signal
        [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
        ev = beh(x).time(locs); % Convert peak locations to seconds  
        ev = ev(ev > r(x,1) & ev < r(x,2)); % extract rewarded trials during post-infusion window
        
        for y = 1:length(beh(x).FP)
            sig = beh(x).FP{y}; % signal that will be aligned to event times
            sig = sig - nanmean(sig); % subtract mean of trace to center on zero
            [mat, time] = getSTA(sig, ev, Fs, [winAcc(1), winAcc(end)]);
            [matBase] = getSTA(sig, ev, Fs, [winBase(1), winBase(end)]);
            matBase = nanmean(matBase,1); % average across baseline window
            cannula(z).a_acc{x,y} = mat - matBase; % save into structure
        end
    end
end

%% PLOT average reward responses
if ~isfield(cannula(1),'a_acc')
    error('No field a_acc - must align photometry to reward delivery before proceeding.');
end
fig = figure; fig.Position(3) = 1000;
for z = 1:2
    a_acc = cannula(choice(z)).a_acc; % extract reward-aligned data from structure
    for y = 1:size(a_acc,2)
        a_mat = []; % initialize temporary output matrix
        for x = 1:size(a_acc,1)
            a_mat(:,x) = nanmean(a_acc{x,y},2); % extract reward-aligned data from structure and average across all rewards
        end
        subplot(1,2,y);
        shadederrbar(time, nanmean(a_mat,2), SEM(a_mat,2), clr{choice(z)}); % plot average across trials
        ylabel(sprintf('%s amplitude (dF/F)',cannula(1).s(1).FPnames{y}));
        xlabel(sprintf('time to %s (s)',xlbl));
        title(sprintf('%s (n = %d): %s vs %s',cannula(1).s(1).FPnames{y},size(a_mat,2),cannula(choice(1)).inf,cannula(choice(2)).inf)); axis square
    end
end
movegui(gcf,'center');

%% Find AMPLITUDE and LATENCY of reward response
if ~isfield(cannula(1),'a_acc')
    error('No field a_acc - must align photometry to reward delivery before proceeding.');
end
fig = figure; fig.Position(3) = 1000;
lag = {}; val = {}; % initialize output
for z = 1:2
    a_acc = cannula(choice(z)).a_acc; % extract reward-aligned data from structure
    Fs = cannula(choice(z)).s(1).Fs; % sampling frequency
    for y = 1:size(a_acc,2)
        r = (winAccPk(y,:)/(1000/Fs) + find(time == 0)); 
        for x = 1:size(a_acc,1)
            pull = a_acc{x,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
            pull = pull(r(1):r(2),:); % extract only segment of time specified as window for reward-related trough
            [a,b] = max(pull); % find local MAXIMUM within range
            c = time(b + r(1) - 1); % convert index to seconds
            c(isnan(a)) = nan; 
            val{y}(x,z) = nanmean(a); % save amplitude of minima
            lag{y}(x,z) = nanmean(c).*1000; % save lag at minima, in milliseconds
        end
        subplot(1,2,y); hold on
        plot(lag{y}(:,z), val{y}(:,z), '.', 'MarkerSize', 20, 'Color', clr{choice(z)}); % plot individual data points for each mouse
        errorbar(nanmean(lag{y}(:,z)), nanmean(val{y}(:,z)), ... % plot average across all mice
        SEM(val{y}(:,z),1), SEM(val{y}(:,z),1), ... % error bars vertical for amplitude
        SEM(lag{y}(:,z),1), SEM(lag{y}(:,z),1), ... % error bars horizontal for latency
        '.', 'MarkerSize', 20, 'Color', clr{choice(z)});
        xlabel(sprintf('time to %s (s)',xlbl)); xlim([winAccPk(y,1) winAccPk(y,2)]);
    end
end
y = 1; subplot(1,2,y); % ACh reward response subplot
    xlabel('time to accel (s)'); xlim([winAccPk(y,1) winAccPk(y,2)]);
    ylabel('ACh peak amplitude (%dF/F)'); ylim([0 10]); yticks([0:2:10]);
    [~,p] = ttest(lag{y}(:,1),lag{y}(:,2)); % statistical test: paired t-test
    [~,p(2)] = ttest(val{y}(:,1),val{y}(:,2)); % statistical test: paired t-test
    title(sprintf('%s vs %s (lag %1.2f, val %1.2f)',cannula(choice(1)).inf,cannula(choice(2)).inf,p(1),p(2)));
    axis square
y = 2; subplot(1,2,y); % DA reward response subplot
    xlabel('time to accel (s)'); xlim([winAccPk(y,1) winAccPk(y,2)]);
    ylabel('DA peak amplitude (%dF/F)'); ylim([0 10]); yticks([0:2:10]);
    [~,p] = ttest(lag{y}(:,1),lag{y}(:,2)); % statistical test: paired t-test
    [~,p(2)] = ttest(val{y}(:,1),val{y}(:,2)); % statistical test: paired t-test
    title(sprintf('%s vs %s (lag %1.2f, val %1.2f)',cannula(choice(1)).inf,cannula(choice(2)).inf,p(1),p(2)));
    axis square
movegui(gcf,'center');