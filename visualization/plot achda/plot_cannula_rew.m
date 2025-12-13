loaded = menu('Already loaded cannula into workspace?','yes','no');
switch loaded
    case 2
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat']); % Select .mat files you want to load
        load(fullfile(fPath,fName));
end

%% INPUTS
nAn = 4;
lickWithin = 1; %CHANGE, lick within this window
winRew = [-1 2]; % CHANGE, window for aligning signal to rewarded trials
winAcc = [-1 1]; % CHANGE, window for aligning signal to acceleration
NumStd = 2; % CHANGE, for immobility trough/peak analysis

%%
% fig = figure;
% plm = 2; pln = nAn;
% clr = {'k','r','b','c','m'};
% lbl = {'acsf','iGluR','nAChR','mAChR','d1d2'};

for z = 1:length(cannula)
    beh = cannula(z).s;
    win = cannula(z).win;
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        rew = beh(x).reward./Fs; % reward delivery times, in seconds
        lick = beh(x).lick./Fs; % lick times, in seconds
        [rewYes, ~, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
        %rewYes: indices of rewarded trials (animal licked within specified window)
        ev = rew(rewYes); % delivery times for rewarded trials, in seconds
        ev = find(ev > win(x,1)*60 & ev < win(x,2)*60); % extract rewarded trials during post-infusion window
        
        for y = 1:length(beh(x).FP)
            sig = beh(x).FP{y}; % signal that will be aligned to event times
            sig = sig - nanmean(sig); % subtract mean of trace to center on zero
            [mat, time] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]);
            cannula(z).a_rew{x,y} = mat; % save into structure
            
%             sp(x) = subplot(plm,pln,x+(y*nAn)); hold on
%             sp(x) = subplot(plm,pln,x+((y-1)*nAn)); hold on
%             shadederrbar(time, nanmean(mat,2), SEM(mat,2), clr{z}); % plot average across trials
%             ylabel(sprintf('%s (%dF/F)',beh(x).FPnames{y}));
%             xlabel('latency to reward (s)');
%             title(sprintf('%s',strtok(beh(x).rec,'-')));
        end
    end
end
%movegui(gcf,'center');

%% Find AMPLITUDE and LATENCY of ACh reward response 
lag = []; val = []; % initialize output
y = 1; % analyzing ACh reward response
r = [200 700]; % range for ACh trough, in milliseconds
r = (r/(1000/Fs) + find(time == 0)); 
for z = 1:length(cannula)
    for x = 1:size(cannula(z).a_rew,1)
        pull = cannula(z).a_rew{x,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
% for z = 1:length(compare)
%     for x = 1:size(compare(z).a_rew,1)
%         pull = compare(z).a_rew{x,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
        pull = pull(r(1):r(2),:); % extract only segment of time specified as window for reward-related trough
        [a,b] = min(pull); % find local minima within range
        c = time(b + r(1) - 1); % convert index to seconds
        c(isnan(a)) = nan; 
        val(x,z) = nanmean(a); % save amplitude of minima
        lag(x,z) = nanmean(c).*1000; % save lag at minima, in milliseconds
    end
end

%% PLOT AMPLITUDE/LATENCY of ACh reward response
figure; hold on
clr = {'k','g'};
for z = 1:2
    plot(lag(:,z), val(:,z), '.', 'MarkerSize', 20, 'Color', clr{z});
    errorbar(nanmean(lag(:,z)), nanmean(val(:,z)), ... 
    SEM(val(:,z),1), SEM(val(:,z),1), ... 
    SEM(lag(:,z),1), SEM(lag(:,z),1), '.', 'MarkerSize', 20, 'Color', clr{z});
end
ylabel('ACh trough amp'); ylim([-5 0]); yticks([-8:2:0]);
xlabel('time to rew (s)'); xlim([200 700]);
[~,p] = ttest2(lag(:,1),lag(:,2));
[~,p(2)] = ttest2(val(:,1),val(:,2));
title(sprintf('ttest2: lag = %1.2f, val = %1.2f',p(1),p(2)));
axis square
movegui(gcf,'center');