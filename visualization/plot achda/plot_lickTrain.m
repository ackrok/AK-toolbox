behLoaded = menu('beh loaded into workspace already?','yes','no');

switch behLoaded
    case 2
        %% Select .mat files you want to add to summary data structu
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
        beh = extractBeh(fPath, fName);
end
if ~exist('beh','var'); error('No variable called beh exists'); end

%% Plot lick rate
lickWithin = 0.25; %CHANGE, lick within this window

figure; plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm);
for x = 1:length(beh)
    Fs = beh(x).Fs;
    rew = beh(x).reward./Fs; % reward delivery times, in seconds
    lick = beh(x).lick./Fs; % lick times, in seconds
    [rewYes, rewNo, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
    %rewYes: indices of rewarded trials (animal licked within specified window)
    %rewNo: indices of non-rewarded trials
    %lickNew: corrected vector of lick times, in seconds
    bin = 0.1; window = [-2 2];
    peth = getClusterPETH(lickNew, rew(rewYes), bin, window); % align lick to rewarded trials
    cts = peth.cts{1}./bin; % licks aligned to rewarded trials, adjust by bin size
    sp(x) = subplot(plm,pln,x); hold on
    plot([0 0],[0 10],'k');
    shadederrbar(peth.time, nanmean(cts,2), SEM(cts,2), 'm'); % plot lick rate
    xlabel('latency to reward (s)');
    ylabel('lick rate (licks/s)'); ylim([0 11])
    title(sprintf('%s: %d%s of %d',beh(x).rec,round(100*length(find(rewYes))/length(rew)),'%',length(rew)));
    %title of subplot is 
    % <recording name (#rewarded trials / total #trials)>
end