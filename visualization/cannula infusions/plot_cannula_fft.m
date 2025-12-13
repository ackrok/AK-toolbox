%% Description
% The code below run FFT analysis and plot output for cannula infusion
% experiments
%
% INPUTS
%   none
%
% OUTPUTS
%   'sub_comp' - cell array with FFT output for all recordings, 
%       each cell contains data for one infusion, subtracting
%       green or red fluorophore, as specified
%   'auc' - area under the curve in 0.5-4Hz frequency band
%       each column contains data for one infusion
%
% Anya Krok, September 2022

%% INPUTS
% winInf = [20 40]; winInf = winInf.*60;
% winInf = [200 800];
winInf = [30 50]; winInf = winInf.*60;

%% LOAD RAW SIGNALS INTO WORKSPACE
loaded = menu('Already loaded raw data into workspace?','yes','yes but update','no');
switch loaded
    case 1
        if ~exist('rawS','var') && ~exist('norm','var')
            [norm, f, flog, p1_mat, rawS] = getCannula_fft(winInf);
        end
    case 2
        if ~exist('rawS','var'); error('No variable called rawS exists'); end
        [norm, f, flog, p1_mat, rawS] = getCannula_fft(winInf, rawS);
    case 3
        [norm, f, flog, p1_mat, rawS] = getCannula_fft(winInf);
end

%% Subtract stable fluorophore (GFP/tdTomato) signal
if ~exist('norm_gfp','var') && ~exist('norm_antd1d2','var')
    [fName_flu,fPath_flu] = uigetfile('*.mat','Select GFP+tdTomato for FFT file');
    load(fullfile(fPath_flu,fName_flu));
end

sub_1 = [norm]; % SUBTRACT
switch rawS(1).fp_lbl
    case 'ACh'; sub_2 = norm_gfp; % FFT ouput: GFP fluorescence signal, average over n = 3 mice
    % case 'DA'; sub_2 = norm_tdt; % FFT ouput: tdTomato fluorescence signal, average over n = 3 mice
    case 'DA'; sub_2 = norm_daAnt; % FFT ouput: rDA fluorescence signal during infusion of DA receptor antagonist, average over n = 4 mice
end
sub_mat = []; for x = 1:size(sub_1,2); sub_mat(:,x) = sub_1(:,x) - nanmean(sub_2,2); end % Subtract avg FFT for mAChR antagonist

%% AUC in specified frequency band
range_auc = [0.5 4]; % Range for calculation of area under the curve, in Hz
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_auc = [find(f_sub == range_auc(1)):find(f_sub == range_auc(2))]; % extract AUC from [0.5 4] Hz
for x = 1:size(sub_mat,2)
    auc(x) = trapz(sub_mat(r_auc,x))/length(r_auc);
end

%% PLOT
ds = 50; % Downsample when plotting so figure is smaller in size
nAn = length(rawS); % Number of recordings and/or animals
fig = figure; fig.Position(3) = 1375;
subplot(1,3,1); hold on
    plot(flog(1:ds:end), sub_mat((1:ds:end),:)); % Plot FFT for each recording
    plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
    plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
    plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
    xlabel('Frequency'); ylabel('Power (a.u.)');
    xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
    legend({rawS.rec});
    title(sprintf('%s FFT %s || win [%d %d]s',rawS(1).fp_lbl,rawS(1).behState,winInf(1),winInf(2))); axis square
subplot(1,3,2); hold on
    shadederrbar(flog(1:ds:end), nanmean(sub_mat((1:ds:end),:),2), SEM(sub_mat((1:ds:end),:),2), 'r'); % Plot FFT averaged across all recordings
    plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
    plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
    plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
    xlabel('Frequency'); ylabel('Power (a.u.)');
    xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
    title(sprintf('%s FFT %s (n = %d)',rawS(1).fp_lbl,rawS(1).behState,nAn)); axis square
subplot(1,3,3); hold on
    j1 = 0.8; j2 = 1.2; jit = j1 + (j2-j1).*rand(nAn,1); % jitter
    plot(jit, auc, '.k', 'MarkerSize', 20); % Plot area under curve power data points for each recording
    errorbar(nanmean(auc), SEM(auc,2), '.r', 'MarkerSize', 20); % Plot area under curve power averaged
    xlim([0.5 1.5]); xticks([1]); xticklabels({'AUC'});
    ylabel('Power (a.u.)'); ylim([0 0.5]);
    title(sprintf('AUC %1.1f-%dHz (mu: %1.2f)',range_auc(1),range_auc(2),nanmean(auc))); axis square
movegui(gcf,'center');