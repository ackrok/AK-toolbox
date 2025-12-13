%% Description
% The code below will use plot pre-analyzed FFT signals. User must have
% already run FFT analysis and have normalized FFT ouput in workspace.
%
% INPUTS
%   'norm' - matrix with FFT output for all recordings
%   'f' - frequency domain vector
%   'flog' - log-scale frequency domain vector
%
% OUTPUTS
%   'sub_mat' - matrix with FFT output for all recordings, subtracting
%               green or red fluorophore, as specified
%   'auc' - area under the curve in 0.5-4Hz frequency band
%
% Anya Krok, July 2022

%% INPUTS
nComp = 2; % Number of comparisons to be made
norm_comp = cell(1,nComp);
norm_comp{1} = [norm(:,[1:4])]; % CHANGE
norm_comp{2} = [norm(:,[5:8])]; % CHANGE
lgd = {'pre','lesion'}; % CHANGE
clr = {'k','g'};
nAn = size(norm_comp{2},2);

%% LOAD pre-analyzed signals into workspace
if ~exist('norm_comp','var')
    [fName,fPath] = uigetfile(['*.mat'],'Load file containing FFT output variable, norm_comp','MultiSelect','Off');
end

%% LOAD fluorophore (GFP/tdTomato) signal
if ~exist('norm_gfp','var') && ~exist('norm_antd1d2','var') && ~exist('norm_tdt','var')
    [fName_flu,fPath_flu] = uigetfile(['*.mat'],'Select GFP+tdTomato for FFT file','MultiSelect','Off');
    load(fullfile(fPath_flu,fName_flu));
end
%
fluorophore = menu(sprintf('Fluorophore FFT signal to use for: %s, %s',lgd{1},lgd{2}),'green','red');
%
switch fluorophore
    case 1; sub = norm_gfp; % FFT ouput: GFP fluorescence signal, average over n = 3 mice
    % case 2; sub_2 = norm_tdt; % FFT ouput: tdTomato fluorescence signal, average over n = 3 mice
    case 2; sub = norm_daAnt; % FFT ouput: rDA fluorescence signal during infusion of DA receptor antagonist, average over n = 4 mice
end
%% SUBTRACT stable fluorophore (GFP/tdTomato) signal
sub_comp = cell(1,nComp);
for y = 1:nComp
    for x = 1:size(norm_comp{y},2)
        sub_comp{y}(:,x) = norm_comp{y}(:,x) - nanmean(sub,2); 
    end
end

%% AUC in specified frequency band
range_auc = [0.5 4]; % Range for calculation of area under the curve, in Hz
auc_comp = cell(1,2);
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_auc = [find(f_sub == range_auc(1)):find(f_sub == range_auc(2))]; % extract AUC from [0.5 4] Hz
for y = 1:length(sub_comp)
    for x = 1:size(sub_comp{y},2)
        auc_comp{y}(x) = trapz(sub_comp{y}(r_auc,x))/length(r_auc);
    end
end

%% PLOT
ds = 50; % Downsample when plotting so figure is smaller in size
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
    for y = 1:length(sub_comp)
        plotme = sub_comp{y}((1:ds:end),:);
        shadederrbar(flog(1:ds:end), nanmean(plotme,2), SEM(plotme,2), clr{y}); % Plot FFT averaged across all recordings
    end
    plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
    plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
    plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
    xlabel('Frequency'); ylabel('Power (a.u.)');
    xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
    title(sprintf('FFT (n = %d)',nAn)); axis square
subplot(1,2,2); hold on
    if diff(cellfun(@length,auc_comp)) == 0
        a = [auc_comp{1}; auc_comp{2}];
        plot(a, '--.k', 'MarkerSize', 20); % Plot matched data points
        [~,p] = ttest(a(1,:),a(2,:)); % Paired test 
    else
        a = nan(2,max(cellfun(@length,auc_comp)));
        for y = 1:2; a(y,[1:length(auc_comp{y})]) = auc_comp{y}; end
        plot(a, '.k', 'MarkerSize', 20); % Plot two-sample data points
        [~,p] = ttest2(a(1,:),a(2,:)); % Two-sample test
    end
    errorbar([0.75 2.25],nanmean(a,2), SEM(a,2), '.', 'MarkerSize', 20, 'Color', clr{2});
    xlim([0.5 2.5]); xticks([1 2]); xticklabels(lgd);
    ylabel('Power (a.u.)'); %ylim([0 0.5]); yticks([0:0.1:0.5]);
    title(sprintf('AUC p = %1.2f',p)); axis square
movegui(gcf,'center');