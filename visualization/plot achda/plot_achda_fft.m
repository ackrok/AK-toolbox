%% Description
% The code below will use FFT to analyze the frequency components present
% in photometry signal with particular emphasis on delta frequency band
% found to be enriched in DLS ACh and DA photometry signals
%
% (1) Dialog window will prompt user to specify inputs. Options are in
%   parentheses!
% (2) If raw signals are not already loaded into workspace in a structure 
%   called rawS, user will then be prompted to select the .mat files that 
%   contain photometry data to be analyzed.
% (3) User will be prompted to select .mat file that contains normalized
%   FFT vectors for recordings of GFP and tdTomato.
% (4) FFT will be run on signal extracted during specified behavioral state,
%   then normalized between 0.1 and 100 Hz, and lastly the FFT power for 
%   GFP or tdTomato will be subtracted.
% (5) FFT output will be plotted if specified to do so.
%
% INPUTS
%   Specify in dialog window.
%   Background functioned needed to run this code: extractEventST
%   Data files needed: raw data files (.mat), GFP+tdTomato data file (.mat)
%
% OUTPUTS
%   'flog' - frequency vector (log-scale),  for plotting FFT output
%   'sub_mat' - matrix with FFT output for all recordings
%   'auc' - area under the curve in 0.5-4Hz frequency band
%
% Anya Krok, July 2022

%% INPUTS
range_norm = [0.1 100]; % Range for normalization of FFT output, in Hz
range_auc = [0.5 4]; % Range for calculation of area under the curve, in Hz
% fluorophore = 'red'; % Options: 'green','red' -- green: use GFP, red: use tdTomato for subtraction
plotYes = 'yes'; % Options: 'yes','no' -- yes plot, no plot

% answer = inputdlg({'Fluorophore: (green, red)','Plot FFT output? (yes, no)'},...
%     'Inputs for FFT', [1 40],{fluorophore,plotYes});
% fluorophore = answer{1};
% plotYes = answer{2};

%% LOAD RAW SIGNALS INTO WORKSPACE
loaded = menu('Already loaded raw data into workspace?','yes','yes but update','no');
switch loaded
    case 2
        rawS = extractRaw_fft(rawS);
    case 3
        %% Select .mat files you want to add to summary data structu
        [fName,fPath] = uigetfile('*.mat','MultiSelect','On');
        rawS = extractRaw_fft(fPath,fName);
end
if ~exist('rawS','var'); error('No variable called rawS exists'); end

%% FFT
[p1_mat, f] = getFft(rawS);
fprintf('FFT done! \n');

%% Normalize FFT
tmpNorm = [];
r = [find(f == range_norm(1)):find(f == range_norm(2))]; % Restrict to [0.01 100]
flog = log10(f(r)); % Log-scale frequency vector
for x = 1:size(p1_mat,2)
    a = log10(p1_mat(r,x));
    a_end = nanmean(log10(p1_mat(find(f == 40):find(f == 50),x)));
    vec_norm = (a - a_end)./(a(1) - a_end); 
    % vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100], scaling so range covers [0 1]
    tmpNorm(:,x) = vec_norm;
end
norm = tmpNorm;
fprintf('Normalization done! \n');

%% Subtract stable fluorophore (GFP/tdTomato) signal
if ~exist('norm_gfp','var') && ~exist('norm_antd1d2','var')
    [fName_flu,fPath_flu] = uigetfile(['*.mat'],'Select GFP+tdTomato for FFT file','MultiSelect','Off');
    load(fullfile(fPath_flu,fName_flu));
end

sub_1 = [norm]; % SUBTRACT
switch rawS(1).fp_lbl
    case 'ACh'; sub_2 = norm_gfp; % FFT ouput: GFP fluorescence signal, average over n = 3 mice
    % case 'ACh'; sub_2 = norm_scop(1:length(flog)); % FFT ouput: ACh fluorescence signal during infusion of mAChR receptor antagonist, average over n = 5 mice    
%     case 'DA'; sub_2 = norm_tdt; % FFT ouput: tdTomato fluorescence signal, average over n = 3 mice
    case 'DA'; sub_2 = norm_daAnt; % FFT ouput: rDA fluorescence signal during infusion of DA receptor antagonist, average over n = 4 mice
end
sub_mat = []; for x = 1:size(sub_1,2); sub_mat(:,x) = sub_1(:,x) - nanmean(sub_2,2); end % Subtract avg FFT for mAChR antagonist

%% AUC in specified frequency band
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_auc = [find(f_sub == range_auc(1)):find(f_sub == range_auc(2))]; % extract AUC from [0.5 4] Hz
for x = 1:size(sub_mat,2)
    auc(x) = trapz(sub_mat(r_auc,x))/length(r_auc);
end

%% PLOT
switch plotYes
    case 'yes'
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
            title(sprintf('%s FFT %s',rawS(1).fp_lbl,rawS(1).behState)); axis square
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
end