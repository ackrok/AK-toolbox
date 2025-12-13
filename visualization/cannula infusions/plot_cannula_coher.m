choice = menu('Cannula structure loaded into workspace?',...
    'YES','No, but will select from pre-saved data file','No, will load from raw data');
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

%% Variables
[plotOpt,ind] = listdlg('PromptString',{'Select All Plots to Generate',...
    'For Multiple: Hold Ctrl and Select'},'ListString',...
    {'DA/ACh cross-corr average traces',...
    'DA/ACh cross-corr statistics',...
    'DA/ACh coherence average traces',...
    'DA/ACh coherence statistics',...
    'DA/ACh phase average traces',...
    'DA/ACh phase statistics'});

%% Cross-correlation Analysis
if any(plotOpt == 1 | plotOpt == 2)
    if ~isfield(cannula,'corr') % If correlation analyis already complete, do not re-run
        for y = 1:length(cannula)
            %%
            beh = cannula(y).s;
            r = cannula(y).win(1,:).*60; % Post infusion window for analysis, in samples
            Fs = beh(1).Fs;
            % Before running analyses, need to overwrite event times for
            % immobility, locomotion, reward to be restricted to post-infusion
            % window for analysis
            [corr_achda, lags, shuff_achda] = AK_corrFP(beh, r(1,:));
            nAn = length(beh); % number of unique animals
            nStates = length(corr_achda); % number of behavioral states
            min_val = nan(nAn,nStates); % initialize output
            min_lag = nan(nAn,nStates); % initialize output
            for z = 1:nStates % iterate over behavioral states
                [min_val(:,z), ii] = min(corr_achda{z,1}(find(lags/Fs == -0.5):find(lags/Fs == 0),:)); % cross-correlation maximal value for each animal during each behavioral state
                ii = ii + find(lags/Fs == -0.5) - 1;
                min_lag(:,z) = lags(ii)./Fs; % latency at minimum, in seconds
                min_lag(isnan(min_val)) = nan;
            end
            cannula(y).corr = corr_achda;
            cannula(y).lags = lags(:)./Fs;
            cannula(y).corrShuff = shuff_achda;
            cannula(y).corrMinVal = min_val;
            cannula(y).corrMinLag = min_lag;
        end
    end
end

%% Coherence/Phase Analysis
if any(plotOpt == 3 | plotOpt == 4 | plotOpt == 5 | plotOpt == 6)
    if ~isfield(cannula,'coher') % If coherence analyis already complete, do not re-run
        for y = 1:length(cannula)
            %%
            beh = cannula(y).s;
            win = cannula(y).win(1,:).*60; % Post infusion window for analysis, in samples
            Fs = beh(1).Fs;
            % Before running analyses, need to overwrite event times for
            % immobility, locomotion, reward to be restricted to post-infusion
            % window for analysis
            [coher_achda, phase_achda, t, f, coher_shuff, phase_shuff] = AK_coherFP(beh, win(1,:));
            nAn = length(beh); % number of unique animals
            nStates = length(coher_achda); % number of behavioral states

            r = [6:42]; % [~,r2] = min(abs(f - 2)); % range for 0.5-4Hz
            coher_avg = nan(nAn,nStates); % initialize output
            phase_avg = nan(nAn,nStates); % initialize output
            for z = 1:nStates % iterate over behavioral states
                coher_avg(:,z) = median(coher_achda{z}(r,:));
                phase_avg(:,z) = median(phase_achda{z}(r,:)); % median within frequency band
            end
            coher_avg(coher_avg == 0) = nan;
            phase_avg(phase_avg == 0) = nan;
            phase_avg = rad2deg(phase_avg); % convert to degrees from radians

            cannula(y).coher = coher_achda;
            cannula(y).phase = phase_achda;
            cannula(y).f = f;
            cannula(y).coherMid = coher_avg;
            cannula(y).phaseMid = phase_avg;
            cannula(y).coherShuff = coher_shuff;
            cannula(y).phaseShuff = phase_shuff;
        end
    end
end

%% Plotting
if ind == 0
    msgbox('Plotting Aborted');
else
    for plotOptNum = 1:length(plotOpt)
        choice = plotOpt(plotOptNum);
        switch choice
            case 1           
                %% PLOT coherence average traces
                fig = figure; hold on; fig.Position([3 4]) = [1375 750];
                clr = {'r','g','b'};
                lbl = {'immobility','locomotion','reward'};
                lags = cannula(1).lags;
                for z = 1:nStates
                    subplot(2,3,z); hold on
                        y = 1; shadederrbar(lags, nanmean(cannula(y).corr{z},2), SEM(cannula(y).corr{z},2), 'k');
                        y = 2; shadederrbar(lags, nanmean(cannula(y).corr{z},2), SEM(cannula(y).corr{z},2), clr{z});
                        plot([0 0],[-1 0.5],'--k');
                        % shadederrbar(lags, nanmean(cannula(1).shuff{z},2), SEM(cannula(1).shuff{z},2), 'k');
                        xlabel('Lag, DA ref (s)'); xlim(win);
                        ylabel('Pearsons r'); ylim([-1 0.5]); yticks([-1:0.2:1]);
                        legend({cannula(1).inf,'',cannula(2).inf},'Location','Southeast');
                        title(sprintf('Correlation - %s',lbl{z})); axis square

                    subplot(2,3,z+3);
                    plotme = [cannula(1).corr{z}, cannula(2).corr{z}]; % concatenate coherence across multiple infusions
                    plotme = plotme(find(lags == win(1)):find(lags == win(2)),:);
                    h = imagesc(lags(find(lags == win(1)):find(lags == win(2))), [1:2*nAn], plotme', [-1 1]); % heatmap of single trials
                    colorbar('off'); colormap(jet(256)); % set color to be jet range
                    xlabel('Lag, DA ref (s)');
                    ylabel('mouse');
                    title(sprintf('top %d: %s | bot %d: %s',nAn,cannula(1).inf,nAn,cannula(2).inf)); axis square
                end
                movegui(gcf,'center');
                
            case 2
                %% DA/ACh cross-corr statistics
                fig = figure; hold on; fig.Position(3) = 1000;
                clr = {'r','g','b'};
                p = [];
                nComp = length(cannula);
                for z = 1:nStates
                    plotme = [cannula(1).corrMinVal(:,z), cannula(2).corrMinVal(:,z)];
                    hold on
                    plot(((nComp*z-1)+[0.25;0.75]).*ones(nComp,nAn), plotme', '.-', 'MarkerSize', 20, 'Color', clr{z});
                    errorbar([z*nComp-1 z*nComp], nanmean(plotme), SEM(plotme,1), '.k', 'MarkerSize', 20);
                    [~,p(z)] = ttest(plotme(:,1),plotme(:,2));
                end
                ylabel('Pearsons r'); ylim([-1 0]); yticks([-1:0.2:0]);
                xlim([0.5 0.5+nComp*nStates]); xticks([1:nComp*nStates]);
                xticklabels(repmat({cannula.inf},1,nStates));
                title(sprintf('Correlation maximum || ttest: imm (%1.3f), loc (%1.3f), rew (%1.3f)',p(1),p(2),p(3)));
                movegui(gcf,'center');
                
            case 3           
                %% PLOT coherence average traces
                fig = figure; hold on; fig.Position([3 4]) = [1375 750];
                clr = {'r','g','b'};
                lbl = {'immobility','locomotion','reward'};
                for z = 1:nStates
                    subplot(2,3,z); hold on
                        y = 1; shadederrbar(f, nanmean(cannula(y).coher{z},2), SEM(cannula(y).coher{z},2), 'k'); % saline
                        y = 2; shadederrbar(f, nanmean(cannula(y).coher{z},2), SEM(cannula(y).coher{z},2), clr{z}); % infusion
                        plot([0.5 0.5],[0 1],'--k'); plot([4 4],[0 1],'--k'); % lines at 0.5, 4Hz
                        % shadederrbar(f, nanmean(cannula(1).coherShuff{z},2), SEM(cannula(1).coherShuff{z},2), 'k');
                        xlabel('Frequency (Hz)');
                        ylabel('Coherence magnitude'); ylim([0 1]); yticks([0:0.2:1]);
                        legend({cannula(1).inf,'',cannula(2).inf});
                        title(sprintf('Coherence - %s',lbl{z})); axis square

                    subplot(2,3,z+3);
                    plotme = [cannula(1).coher{z}, cannula(2).coher{z}]; % concatenate coherence across multiple infusions
                    h = imagesc(f, [1:2*nAn], plotme', [0 1]); % heatmap of single trials
                    colorbar('off'); colormap(jet(256)); % set color to be jet range
                    xlabel('Frequency (Hz)');
                    ylabel('mouse');
                    title(sprintf('top %d: %s | bot %d: %s',nAn,cannula(1).inf,nAn,cannula(2).inf)); axis square
                end
                movegui(gcf,'center');

            case 4
                %% PLOT coherence statistics
                fig = figure; hold on; fig.Position(3) = 1000;
                clr = {'r','g','b'};
                p = [];
                nComp = length(cannula);
                for z = 1:nStates
                    plotme = [cannula(1).coherMid(:,z), cannula(2).coherMid(:,z)];
                    hold on
                    plot(((nComp*z-1)+[0.25;0.75]).*ones(nComp,nAn), plotme', '.-', 'MarkerSize', 20, 'Color', clr{z});
                    errorbar([z*nComp-1 z*nComp], nanmean(plotme), SEM(plotme,1), '.k', 'MarkerSize', 20);
                    [~,p(z)] = ttest(plotme(:,1),plotme(:,2));
                end
                ylabel('Coherence magnitude'); ylim([0 1]); yticks([0:0.2:1]);
                xlim([0.5 0.5+nComp*nStates]); xticks([1:nComp*nStates]);
                xticklabels(repmat({cannula.inf},1,nStates));
                title(sprintf('Coherence maximum || ttest: imm (%1.3f), loc (%1.3f), rew (%1.3f)',p(1),p(2),p(3)));
                movegui(gcf,'center');
                
            case 5
                %% PLOT phase offset average traces
                fig = figure; hold on; fig.Position([3 4]) = [1375 750];
                clr = {'r','g','b'};
                lbl = {'immobility','locomotion','reward'};
                for z = 1:nStates
                    subplot(2,3,z); hold on
                        y = 1; shadederrbar(f, nanmean(rad2deg(cannula(y).phase{z}),2), SEM(rad2deg(cannula(y).phase{z}),2), 'k'); % saline
                        y = 2; shadederrbar(f, nanmean(rad2deg(cannula(y).phase{z}),2), SEM(rad2deg(cannula(y).phase{z}),2), clr{z}); % infusion
                        plot([0.5 0.5],[0 180],'--k'); plot([4 4],[0 180],'--k'); % lines at 0.5, 4Hz
                        % shadederrbar(f, nanmean(cannula(1).phaseShuff{z},2), SEM(cannula(1).phaseShuff{z},2), 'k');
                        xlabel('Frequency (Hz)');
                        ylabel('Phase offset (deg)'); ylim([0 180]); yticks([0:45:180]);
                        legend({cannula(1).inf,'',cannula(2).inf});
                        title(sprintf('Phase - %s',lbl{z})); axis square

                    subplot(2,3,z+3);
                    plotme = [cannula(1).phase{z}, cannula(2).phase{z}]; % concatenate phase offset across multiple infusions
                    plotme = rad2deg(plotme);
                    h = imagesc(f, [1:2*nAn], plotme', [-180 180]); % heatmap of single trials
                    colorbar('off'); colormap(jet(256)); % set color to be jet range
                    xlabel('Frequency (Hz)');
                    ylabel('mouse');
                    title(sprintf('top %d: %s | bot %d: %s',nAn,cannula(1).inf,nAn,cannula(2).inf)); axis square
                end
                movegui(gcf,'center');
                
            case 6
                %% PLOT coherence statistics
                fig = figure; hold on; fig.Position(3) = 1000;
                clr = {'r','g','b'};
                p = [];
                nComp = length(cannula);
                for z = 1:nStates
                    plotme = [cannula(1).phaseMid(:,z), cannula(2).phaseMid(:,z)];
                    hold on
                    plot(((nComp*z-1)+[0.25;0.75]).*ones(nComp,nAn), plotme', '.-', 'MarkerSize', 20, 'Color', clr{z});
                    errorbar([z*nComp-1 z*nComp], nanmean(plotme), SEM(plotme,1), '.k', 'MarkerSize', 20);
                    [~,p(z)] = ttest(plotme(:,1),plotme(:,2));
                end
                ylabel('Phase offset (deg)'); ylim([0 180]); yticks([0:45:180]);
                xlim([0.5 0.5+nComp*nStates]); xticks([1:nComp*nStates]);
                xticklabels(repmat({cannula.inf},1,nStates));
                title(sprintf('Phase Offset || ttest: imm (%1.3f), loc (%1.3f), rew (%1.3f)',p(1),p(2),p(3)));
                movegui(gcf,'center');
        end
    end
end