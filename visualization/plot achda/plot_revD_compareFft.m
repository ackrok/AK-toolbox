% Description
% The code below will use compare frequency bands enriched in specified
% photometry signals during specified behavioral state
%
% INPUTS
%   'norm' - matrix with FFT output for all recordings
%   'f' - frequency domain vector
%   'flog' - log-scale frequency domain vector
%
% OUTPUTS
%   'auc' - area under the curve in 0.5-4Hz frequency band
%
% Anya Krok, July 2022
%
%% DIRECTIONS
% If raw photometry signals NOT loaded into workspace,
% (1) 1st pop-up window: select PRE-lesion data files
% (2) 2nd pop-up window: select LESION data files for same mice
% (3) select gfp, tdtomato files if not already loaded into workspace
% Need to select ones for corresponding mice, but the code will ensure that
% all pairs have same mouse ID

%% INPUTS
choice = menu('Raw photometry signal loaded into workspace?',...
    'yes','no');
switch choice
    case 2
        raw1 = extractRaw_fft; for x = 1:length(raw1); raw1(x).inf = 'pre'; end
        raw2 = extractRaw_fft; for x = 1:length(raw2); raw2(x).inf = 'lesion'; end
        rawS = [raw1, raw2];
        [p1_mat, f] = getFft(rawS);
        tmpNorm = [];
        range_norm = [0.1 100]; % Range for normalization of FFT output, in Hz
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
end

%% EXTRACT DATA for comparison
switch choice
    case 2
        rec = {}; for x = 1:length(rawS); rec{x} = strtok(rawS(x).rec,'-'); end % Extract mouse names
        uniAn = unique(rec); % Unique mouse
        uniInf = unique({rawS.inf}); uniInf = fliplr(uniInf); % Unique infusions
        nComp = length(uniInf); % Number of comparisons to be made
        norm_comp = cell(1,nComp);
        idxStore = cell(1,nComp);
        for x = 1:nComp
            idx = find(strcmp({rawS.inf}, uniInf{x}));
            norm_comp{x} = norm(:,idx); % Extract normalized FFT ouput for matching infusions into single cell array
            idxStore{x} = idx;
        end
        lbl_inf = uniInf; % Store unique infusion labels
        lbl_fp = rawS(1).fp_lbl; % Store photometry label
        nAn = length(uniAn); % Number of unique animals
        for x = 1:nAn % Check that mouse ID is the same across data in different infusion conditions before proceeding
            tmpID = cell(1,nComp); for y = 1:nComp; tmpID{y} = rec{idxStore{2}(x)}; end
            if ~all(strcmp(tmpID, rec{idxStore{1}(x)})) % If all mouse ID do not match ID of this data in first infusion condition
                error('Data across infusions does not have matching mouse IDs. \n Re-load data, making sure to select data files with matching mouse IDs. \n',x);
            end
        end
end

%% LOAD pre-analyzed signals into workspace
if ~exist('norm_comp','var')
    [fName,fPath] = uigetfile('*.mat', 'Select data file that contains norm_comp f flog lbl_inf lbl_fp');
    if ~exist('norm_comp','var')
        error('No variable norm_comp in workspace.');
    end
end

%% LOAD fluorophore (GFP/tdTomato) signal
if ~exist('norm_gfp','var') && ~exist('norm_antd1d2','var') && ~exist('norm_tdt','var')
    [fName_flu,fPath_flu] = uigetfile('*.mat','Select GFP+tdTomato for FFT file');
    load(fullfile(fPath_flu,fName_flu));
end
% fluorophore = menu('Fluorophore FFT signal to subtract','green','red');
switch lbl_fp
    case 'ACh'; sub = norm_gfp; % FFT ouput: GFP fluorescence signal, average over n = 3 mice
    % case 'ACh'; sub = norm_achAnt; % FFT ouput: ACh fluorescence signal during infusion of mAChR receptor antagonist, average over n = 5 mice    
    % case 'DA'; sub = norm_tdt; % FFT ouput: tdTomato fluorescence signal, average over n = 3 mice
    case 'DA'; sub = norm_daAnt; % FFT ouput: rDA fluorescence signal during infusion of DA receptor antagonist, average over n = 4 mice
end

%% SUBTRACT stable fluorophore (GFP/tdTomato) signal
sub_comp = cell(1,length(norm_comp));
for y = 1:length(norm_comp)
    for x = 1:size(norm_comp{y},2)
        sub_comp{y}(:,x) = norm_comp{y}(:,x) - nanmean(sub,2); 
    end
end
[~,b] = cellfun(@size, sub_comp);
nAn = max(b);
nComp = length(sub_comp);

%% AUC in specified frequency band
range_auc = [0.5 4]; % Range for calculation of area under the curve, in Hz
auc = nan(size(sub_comp{1},2), length(sub_comp));
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_auc = [find(f_sub == range_auc(1)):find(f_sub == range_auc(2))]; % extract AUC from [0.5 4] Hz
for y = 1:length(sub_comp)
    for x = 1:size(sub_comp{y},2)
        auc(x,y) = trapz(sub_comp{y}(r_auc,x))/length(r_auc);
    end
end
auc(auc == 0) = nan;
auc(auc < 0) = 0;

%% PLOT
fig = figure; fig.Position(3) = 1000;
ds = 50; % Downsample when plotting so figure is smaller in size
clr = {'k','g'};
subplot(1,2,1); hold on
    for y = 1:nComp
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
    a = auc';
    plot(a, '--.k', 'MarkerSize', 20); % Plot matched data points
    [~,p] = ttest(a(1,:),a(2,:)); % Paired test 
    errorbar([0.75 2.25],nanmean(a,2), SEM(a,2), '.', 'MarkerSize', 20, 'Color', clr{2});
    xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl_inf);
    ylabel('Power (a.u.)'); ylim([0 0.7]); yticks([0:0.1:0.7]);
    title(sprintf('AUC p = %1.2f',p)); axis square
movegui(gcf,'center');