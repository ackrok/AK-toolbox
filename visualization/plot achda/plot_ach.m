behLoaded = menu('beh loaded into workspace already?','yes','no');

switch behLoaded
    case 2
        %% Select .mat files you want to add to summary data structu
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
        if ~iscell(fName); fName = {fName}; end

        %% Extract data
        beh = struct;
        for y = 1:length(fName) 
            load(fullfile(fPath,fName{y})); 
            [an,j2] = strtok(fName{y},'_'); day = strtok(j2,'_');
            x = 1+length(beh);
            beh(x).rec = [an,'-',day]; 
            beh(x).site = 'DLS';
            beh(x).task = 'reward';

            %% Photometry
            beh(x).Fs = data.gen.Fs; % Sampling frequency, in Hz
            beh(x).time = data.final.time; % Time vector
            beh(x).FP = data.final.FP;
            beh(x).nbFP = data.final.nbFP; % Photometry signal(s)
            beh(x).FPnames = data.final.FPnames; % Names of photometry signal(s)

            %% Movement
            if isfield(data.final,'vel') % If movement data exists
                [~,ii] = max(abs(data.final.vel));
                if data.final.vel(ii) < 0
                    data.final.vel = -data.final.vel; % Flip velocity for recordings on IV rig#1, which has inverted positional encoder signal
                    save(fullfile(fPath,fName{z}),'data'); % Overwrite data file to adjust
                end
                beh(x).vel = data.final.vel; % Velocity signal
                beh(x).on = data.final.mov.onsets; beh(x).off = data.final.mov.offsets;                 % Movement onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
                beh(x).onRest = data.final.mov.onsetsRest; beh(x).offRest = data.final.mov.offsetsRest; % Rest onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
            else % Open field
                try 
                    beh(x).cam = data.final.cam.on; 
                catch
                    signal = data.acq.opto.trace;
                    sigEdge = data.gen.params.FP.sigEdge; 
                    rawFs = data.gen.params.acqFs; dsRate = data.gen.params.dsRate;
                    if sigEdge ~= 0
                        signal = signal((sigEdge*rawFs)+1:end-(sigEdge*rawFs));
                    end
                    camOn = getPulseOnsetOffset (signal, 0.5);
                    camOn_Fs = round(camOn/dsRate);
                    beh(x).cam = camOn_Fs;
                end
            end

            %% Lick/Reward
            if isfield(data.acq,'rew')
                % data = processReward(data, data.gen.params);
                beh(x).task = 'reward';
                beh(x).reward = data.final.rew.onset;    % Reward delivery time in sampling freq (data.gen.Fs), NOT in seconds
                beh(x).lick = data.final.lick.onset;        % Lick times in sampling freq (data.gen.Fs), NOT in seconds
                beh(x).lickVec = data.final.lick.trace;        % Lick trace
            end

            %%
            fprintf('Extracted from %s\n',fName{y});
        end
        if isempty(beh(1).Fs); beh(1) = []; end
end

%% Variables
[plotOpt,ind] = listdlg('PromptString',{'Select All Plots to Generate',...
    'For Multiple: Hold Ctrl and Select'},'ListString',...
    {'Signal over time',...
    'Align to Reward',...
    'Align to Reward+Single Trial',...
    'Align to Reward+1st Lick',...
    'Align to Acceleration',...
    'Immobility Troughs/Peaks'});
%
y = 1; %only ACh
plotSingleTrial = 0; % CHANGE, if = 0 then will not plot single trial data;
lickWithin = 0.25; %CHANGE, lick within this window
winRew = [-1 2]; % CHANGE, window for aligning signal to rewarded trials
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
    figure; y = 1;
    for x = 1:length(beh)
        subplot(length(beh),1,x); hold on
        plot(beh(x).time, beh(x).FP{y}, 'k');
        plot(beh(x).time, -5+getAcc(beh(x).vel),'r');
        ylabel('%dF/F')
        title(sprintf('%s',beh(x).rec));
    end
    subplot(length(beh),1,1); legend({beh(x).FPnames{y}, 'acc'});
    subplot(length(beh),1,length(beh)); xlabel('Time (s)'); 
    movegui(gcf,'center');

case {2, 3}
    %% Plot photometry to rewarded trials
    % lickWithin = 0.25; %CHANGE, lick within this window
    % plotSingleTrial = 1; % CHANGE, if = 0 then will not plot single trial data;
    % y = 1; %only ACh
    % winRew = [-1 2]; % CHANGE, window for aligning signal to events

    figure; 
    switch choice % set subplot # for figure
        case 2; plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm); % if not plotting single trials
        case 3; plm = 2; pln = length(beh); % if yes plotting single trials
    end
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        rew = beh(x).reward./Fs; % reward delivery times, in seconds
        lick = beh(x).lick./Fs; % lick times, in seconds
        [rewYes, ~, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
        %rewYes: indices of rewarded trials (animal licked within specified window)
        ev = rew(rewYes); % CHANGE, aliging signal to rewarded trails

        sig = beh(x).FP{y}; % signal that will be aligned to event times
        sig = sig - nanmean(sig); % subtract mean of trace to center on zero
        [mat, time] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]);

        sp(x) = subplot(plm,pln,x); hold on
        plot([0 0],[-4 2],'k'); % plot line at reward delivery t = 0
        shadederrbar(time, nanmean(mat,2), SEM(mat,2), 'g'); % plot average across trials
        xlabel('latency to reward (s)'); ylabel('ACh3.0 (%dF/F)');
        title(sprintf('%s (%d trials)',beh(x).rec,length(find(rewYes))),'Interpreter','none');
        %title of subplot is <recording name (#rewarded trials)>

        if choice == 3
            subplot(plm,pln,x+length(beh)); hold on
            h = imagesc(time, [1:length(ev)], mat', [-5 15]); % heatmap of single trials
            colorbar('eastoutside'); colormap(jet(256)); % set color to be jet range
            if x ~= length(beh); colorbar('hide'); end
            
            title('single trial');
            xlabel('latency (s)'); xlim([-1 2]);
            ylabel('trial #'); ylim([0 length(ev)]);
            h2 = h.Parent; % h2.CLim = [-0.2 0.6];
        end
    end
    linkaxes(sp,'y');
    movegui(gcf,'center');


case 5
    %% Plot photometry to peaks of positive acceleration
    % winAcc = [-1 1];
    figure;
    plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm);
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        acc = getAcc(beh(x).vel); % Extract acceleration signal
        [~,locs] = findpeaks(acc,'MinPeakProminence',0.5,'MinPeakDistance',0.5); % Location of peaks, using findpeaks function
        ev = beh(x).time(locs); % Convert peak locations to seconds  

        sig = beh(x).FP{y}; % signal that will be aligned to event times
        sig = sig - nanmean(sig); % subtract mean of trace to center on zero
        [mat, time] = getSTA(sig, ev, Fs, [winAcc(1), winAcc(end)]);

        sp(x) = subplot(plm,pln,x); hold on
        plot([0 0],[-5 10],'k'); % plot line at reward delivery t = 0
        shadederrbar(time, nanmean(mat,2), SEM(mat,2), 'g'); % plot average across trials
        xlabel('latency to accel (s)'); ylabel('ACh3.0 (%dF/F)');
        title(sprintf('%s',beh(x).rec),'Interpreter','none');
        %title of subplot is <recording name (#rewarded trials)>
    end
    linkaxes(sp,'y');
    movegui(gcf,'center');


case 6
    %% Plot immobility pause peak
    % NumStd = 2;
    [amp, dur, freq, thres] = getImmPausePeak(beh);
    freq(isnan(freq)) = 0;
    dur(isnan(dur)) = 0; dur = dur.*(1000/Fs); % Adjust from samples to ms
    amp(isnan(amp)) = 0; amp = abs(amp); % adjust to be absolute amplitude
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

case 4
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

        sig = beh(x).FP{y}; % signal that will be aligned to event times
        sig = sig - nanmean(sig); % subtract mean of trace to center on zero
        [mat, time] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]);

        sp(x) = subplot(plm,pln,x); hold on
        plot([0 0],[-4 2],'k'); % plot line at reward delivery t = 0
        shadederrbar(time, nanmean(mat,2), SEM(mat,2), 'g'); % plot average across trials
        xlabel('latency (s)'); ylabel('ACh3.0 (%dF/F)');
        title(sprintf('%s (%d trials)',beh(x).rec,length(find(rewYes))),'Interpreter','none');
        %title of subplot is <recording name (#rewarded trials)>
        
        % % aligning photometry to FIRST lick
        bin = 1/1000;
        peth = getClusterPETH(lickNew, ev, bin, [0 1]); % PETH: lick aligned to rewarded trials in 1 ms bins
        cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
        [~, lickFirst] = max(cts~=0, [], 1); % Find first non-zero index for each trial
        lickFirst = lickFirst(:).*bin + ev(:); % First lick after delivery for rewarded trials
        [mat, time] = getSTA(sig, lickFirst, Fs, [winRew(1), winRew(end)]);
        % plot([0 0],[-4 2],'k'); % plot line at reward delivery t = 0
        shadederrbar(time, nanmean(mat,2), SEM(mat,2), 'b'); % plot average across trials
%         xlabel('latency to 1st lick (s)'); ylabel('ACh3.0 (%dF/F)');
%         title('1st lick');
        
        % % aligning photometry to ALL licks
%         subplot(plm,pln,x+length(beh)); hold on
%         [mat, time] = getSTA(sig, lickNew, Fs, [winRew(1), winRew(end)]);
%         shadederrbar(time, nanmean(mat,2), SEM(mat,2), 'g'); % plot average across trials
%         xlabel('latency to lick (s)'); ylabel('ACh3.0 (%dF/F)');
%         title('licks');
    end
    linkaxes(sp,'y');
    %%
    end
    end
end