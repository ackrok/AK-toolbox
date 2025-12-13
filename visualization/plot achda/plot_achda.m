behLoaded = menu('beh loaded into workspace already?','yes','no');

switch behLoaded
    case 2
        %% Select .mat files you want to add to summary data structu
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
        beh = extractBeh(fPath, fName);
end
if ~exist('beh','var'); error('No variable called beh exists'); end

%% Variables
[plotOpt,ind] = listdlg('PromptString',{'Select All Plots to Generate',...
    'For Multiple: Hold Ctrl and Select'},'ListString',...
    {'Signal over time',...
    'Non-baselined signal',...
    'Align to Reward',...
    'Align to Reward+Single Trial',...
    'Align to Reward+1st Lick',...
    'Align to Acceleration',...
    'Immobility Troughs/Peaks'});
%
plotSingleTrial = 0; % CHANGE, if = 0 then will not plot single trial data;
lickWithin = 0.25; %CHANGE, lick within this window
winRew = [-1 2]; % CHANGE, window for aligning signal to rewarded trials
winBase = [-4 -1]; % window for baseline
winAcc = [-1 1]; % CHANGE, window for aligning signal to acceleration
NumStd = 2; % CHANGE, for immobility trough/peak analysis

% answer = inputdlg({...
%     'Plot single trial data? (1 is yes, 0 is no)',...
%     'Rewarded trials = lick within Xs',...
%     'numSTD for trough/peak analysis'},...
%     'Input', 1,...
%     {num2str(plotSingleTrial),num2str(lickWithin),num2str(NumStd)});
% plotSingleTrial = str2num(answer{1});
% lickWithin = str2num(answer{2});
% NumStd = str2num(answer{3});

%%
if ind == 0
    msgbox('Plotting Aborted');
else
    for plotOptNum = 1:length(plotOpt)
        choice = plotOpt(plotOptNum);
        switch choice

case 1
    %% Plot photometry
    figure; 
    for x = 1:length(beh)
        subplot(length(beh),1,x); hold on
        Fs = beh(x).Fs;
        for y = 1:length(beh(x).FP); clr = {'g','m'};
            plot(beh(x).time, beh(x).FP{y}, clr{y});
        end
        if isfield(beh,'vel')
            plot(beh(x).time, -5+getAcc(beh(x).vel),'r');
        end
        if isfield(beh,'reward')
            rew = beh(x).reward./Fs; % reward delivery times, in seconds
            lick = beh(x).lick./Fs; % lick times, in seconds
            plot([rew,rew]', [0;10].*ones(2,length(rew)), 'b', 'LineWidth', 1.5);
            [~, ~, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
            plot([lickNew,lickNew]', [-4;-2].*ones(2,length(lickNew)), 'k');
        end
        ylabel('%dF/F')
        title(sprintf('%s',beh(x).rec));
    end
    subplot(length(beh),1,1); legend(beh(x).FPnames);
    subplot(length(beh),1,length(beh)); xlabel('Time (s)'); 
    movegui(gcf,'center');
    
case 2
    %% Plot non-baselined photometry
    clr = {'g','m'};
    fig = figure; fig.Position(3) = 1000;
    for x = 1:length(beh)
        for y = 1:length(beh(x).nbFP)
            subplot(length(beh), 2, 2*x + (y-2));
            plot(beh(x).time, beh(x).nbFP{y}, clr{y});
            title(sprintf('%s n.b. %s', beh(x).rec, beh(x).FPnames{y}));
        end
    end
    
case {3, 4}
    %% Plot photometry to rewarded trials
    % lickWithin = 0.25; %CHANGE, lick within this window
    % plotSingleTrial = 1; % CHANGE, if = 0 then will not plot single trial data;
    % y = 1; %only ACh
    % winRew = [-1 2]; % CHANGE, window for aligning signal to events

    figure; 
    switch choice % set subplot # for figure
        case 3; plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm); % if not plotting single trials
        case 4; plm = 3; pln = length(beh); % if yes plotting single trials
    end
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        rew = beh(x).reward./Fs; % reward delivery times, in seconds
        lick = beh(x).lick./Fs; % lick times, in seconds
        [rewYes, ~, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
        %rewYes: indices of rewarded trials (animal licked within specified window)
        ev = rew(rewYes); % CHANGE, aliging signal to rewarded trails
        
        sp(x) = subplot(plm,pln,x); hold on
            plot([0 0],[-4 2],'Color',[0 0 0 0.2]); % plot line at reward delivery t = 0
            
        for y = 1:length(beh(x).FP); clr = {'g','m'};
            sig = beh(x).FP{y}; % signal that will be aligned to event times
            sig = sig - nanmean(sig); % subtract mean of trace to center on zero
            [mat, time] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]);
            [matBase] = getSTA(sig, ev, Fs, [winBase(1), winBase(end)]);
            matBase = nanmean(matBase,1); % average across baseline window
            shadederrbar(time, nanmean(mat - matBase,2), SEM(mat - matBase,2), clr{y}); % plot average across trials
        end
        xlabel('latency to reward (s)'); ylabel('FP (%dF/F)');
        title(sprintf('%s (%d trials)',beh(x).rec,length(find(rewYes))),'Interpreter','none');
        %title of subplot is <recording name (#rewarded trials)>

        if choice == 4
            for y = 1:length(beh(x).FP)
                subplot(plm,pln,x+(y*length(beh))); hold on
                sig = beh(x).FP{y}; % signal that will be aligned to event times
                sig = sig - nanmean(sig); % subtract mean of trace to center on zero
                [mat, time] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]);
                h = imagesc(time, [1:length(ev)], mat', [-5 15]); % heatmap of single trials
                colorbar('eastoutside'); colormap(jet(256)); % set color to be jet range
                if x ~= length(beh); colorbar('hide'); end
            
                title(sprintf('%s',beh(x).FPnames{y}));
                xlabel('latency (s)'); xlim([-1 2]);
                ylabel('trial #'); ylim([0 length(ev)]);
                h2 = h.Parent; % h2.CLim = [-0.2 0.6];
            end
        end
    end
    linkaxes(sp,'y');
    movegui(gcf,'center');
    
case 5
    %% Photometry to Lick
    % lickWithin = 0.25; %CHANGE, lick within this window
    % plotSingleTrial = 1; % CHANGE, if = 0 then will not plot single trial data;
    % y = 1; %only ACh
    % winRew = [-1 2]; % CHANGE, window for aligning signal to events

    figure; 
    plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm);
    % plm = 2; pln = length(beh);
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        rew = beh(x).reward./Fs; % reward delivery times, in seconds
        lick = beh(x).lick./Fs; % lick times, in seconds
        [rewYes, ~, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
        %rewYes: indices of rewarded trials (animal licked within specified window)
        ev = rew(rewYes); % CHANGE, aliging signal to rewarded trails
        
        % % aligning photometry to FIRST lick
        bin = 1/1000;
        peth = getClusterPETH(lickNew, ev, bin, [0 1]); % PETH: lick aligned to rewarded trials in 1 ms bins
        cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
        [~, lickFirst] = max(cts~=0, [], 1); % Find first non-zero index for each trial
        lickFirst = lickFirst(:).*bin + ev(:); % First lick after delivery for rewarded trials
        
        sp(x) = subplot(plm,pln,x); hold on
            plot([0 0],[-4 2],'Color',[0 0 0 0.2]); % plot line at reward delivery t = 0
            
        for y = 1:length(beh(x).FP); clr = {'g','m'}; clr2 = {'b','r'};  
            sig = beh(x).FP{y}; % signal that will be aligned to event times
            sig = sig - nanmean(sig); % subtract mean of trace to center on zero
            [mat, time] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]); % aligning photometry to reward DELIVERY
            shadederrbar(time, nanmean(mat,2), SEM(mat,2), clr{y}); % plot average across trials
            [mat, time] = getSTA(sig, lickFirst, Fs, [winRew(1), winRew(end)]); % aligning photometry to FIRST lick
            shadederrbar(time, nanmean(mat,2), SEM(mat,2), clr2{y}); % plot average across trials
        end
        xlabel('latency (s)'); ylabel('FP (%dF/F)');
        title(sprintf('%s (%d trials)',beh(x).rec,length(find(rewYes))),'Interpreter','none');
        %title of subplot is <recording name (#rewarded trials)>
        
        % % aligning photometry to ALL licks
%         subplot(plm,pln,x+length(beh)); hold on
%         [mat, time] = getSTA(sig, lickNew, Fs, [winRew(1), winRew(end)]);
%         shadederrbar(time, nanmean(mat,2), SEM(mat,2), 'g'); % plot average across trials
%         xlabel('latency to lick (s)'); ylabel('ACh3.0 (%dF/F)');
%         title('licks');
    end
    linkaxes(sp,'y');

case 6
    %% Plot photometry to peaks of positive acceleration
    % winAcc = [-1 1];
    figure;
    plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm);
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        acc = getAcc(beh(x).vel); % Extract acceleration signal
        [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
        ev = beh(x).time(locs); % Convert peak locations to seconds  
        
        sp(x) = subplot(plm,pln,x); hold on
            plot([0 0],[-1 5],'Color',[0 0 0 0.2]); % plot line at reward delivery t = 0
            
        for y = 1:length(beh(x).FP); clr = {'g','m'};
            sig = beh(x).FP{y}; % signal that will be aligned to event times
            sig = sig - nanmean(sig); % subtract mean of trace to center on zero
            [mat, time] = getSTA(sig, ev, Fs, [winAcc(1), winAcc(end)]);
            [matBase] = getSTA(sig, ev, Fs, [winBase(1), winBase(end)]);
            matBase = nanmean(matBase,1); % average across baseline window
            shadederrbar(time, nanmean(mat - matBase,2), SEM(mat - matBase,2), clr{y}); % plot average across trials
        end
        xlabel('latency to accel (s)'); ylabel('FP (%dF/F)');
        title(sprintf('%s',beh(x).rec),'Interpreter','none');
        %title of subplot is <recording name (#rewarded trials)>
    end
    linkaxes(sp,'y');
    movegui(gcf,'center');

case 7
    %% Plot immobility pause peak
    % NumStd = 2;
    [amp, dur, freq, thres] = getImmPausePeak(beh);
    Fs = beh(1).Fs;
    freq(isnan(freq)) = nan;
    dur(isnan(dur)) = nan; dur = dur.*(1000/Fs); % Adjust from samples to ms
    amp(isnan(amp)) = nan; amp = abs(amp); % adjust to be absolute amplitude
    plotme = cell(1,3); plotme{1} = freq; plotme{2} = dur; plotme{3} = amp;
    
    fig = figure; fig.Position(3) = 1375;
    lbl = {'trough','peak'};
    nAn = length(beh);
    for ii = 1:length(plotme)
        a = plotme{ii};
        subplot(1,3,ii); hold on
        jit = []; 
        j1 = 0.8; j2 = 1.2; jit(:,1) = j1 + (j2-j1).*rand(nAn,1); % jitter for pause
        j1 = 1.8; j2 = 2.2; jit(:,2) = j1 + (j2-j1).*rand(nAn,1); % jitter for peak
        plot(jit,a,'.k','MarkerSize',20); % plot raw values per animal
        errorbar([1,2],nanmean(a),SEM(a,1),'.r','MarkerSize',20); % plot average with error bars across animals
        xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl);
        axis square;
    end

    ii = 1; subplot(1,3,ii); % frequency of pauses and peaks
        a = plotme{ii};
        ylabel('frequency (events/s)'); ylim([0 1]); yticks([0:0.5:1]);
        title(sprintf('imm freq (trough mu: %1.2f | peak mu: %1.2f)',(nanmean(a(:,1))),(nanmean(a(:,2))))); % title includes average for peaks and pause
    ii = 2; subplot(1,3,ii); % duration of pause and peaks
        a = plotme{ii};    
        ylabel('duration (ms)'); ylim([0 400]); yticks([0:100:500]);
        title(sprintf('duration (trough mu: %d | peak mu: %d)',round(nanmean(a(:,1))),round(nanmean(a(:,2))))); % title includes average for peaks and pause
    ii = 3; subplot(1,3,ii); % amplitude of pause and peaks
        a = plotme{ii};
        ylabel('abs amplitude (%dF/F)'); ylim([0 6]); yticks([0:2:6]);
        title(sprintf('amplitude (trough mu: %1.1f | peak mu: %1.1f)',(nanmean(a(:,1))),(nanmean(a(:,2))))); % title includes average for peaks and pause

    movegui(gcf,'center');

    %%
    end
    end
end