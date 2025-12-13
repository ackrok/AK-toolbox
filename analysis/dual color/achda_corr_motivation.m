%Cross-correlation between ACh and DA photometry signals acquired during
%dual recordings
%
%   Description: This script is for running cross-correlation analysis on
%   two continuous photometry signals (ACh, DA) after isolating time
%   periods when animal is immobile (no movement, no reward)
%
%   Author: Anya Krok, December 2021

%% 
% good = [1:length(modAChDA)]; rmv = [3 5 18 22:24 27:28 30 33:34 36 41:42]; good(rmv) = []; % acceleration
% good = [10:16,20,22:29,32:35,40:41,44:47]; % reward
good_rew = [10:16,20,22:26,28:29,32:35,40:42,44:46];
beh = modAChDA(good_rew);

% beh = dual;

%%
mat = struct;
corr_cell = cell(3,2);
for a = 1:3; for b = 1:2; corr_cell{a,b} = nan(501,length(beh)); end; end
corr_shuff = cell(3,3);

h = waitbar(0, 'cross corr');
for x = 1:length(beh); y = [2 1]; %CHANGE - which FP signal to run CCG or ACG on
    
    %% extract signals
    fp_mat = [];
    fp_mat(:,1) = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure
    fp_mat(:,1) = fp_mat(:,1) - nanmean(fp_mat(:,1)); % subtract baseline (mean of entire photometry signal) from fp
    fp_mat(:,2) = beh(x).FP{y(2)}; % dopamine
    fp_mat(:,2) = fp_mat(:,2) - nanmean(fp_mat(:,2)); % subtract baseline (mean of entire photometry signal) from fp
    % fp_mat(:,2) = filterFP(fp_mat(:,2),Fs,[0.5 4],10,'bandpass'); 
    
    %% Behavioral states
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+50, 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    
    %% Divide data to assess early/late motivational state
    idx_cell = cell(3,2);
    segment = 3; % how many segments to divide recording into
    segment = reshape([1:length(fp_mat)],[length(fp_mat)/segment, segment]);
    segment(:,2) = []; % remove center segment, keep only early and late
    for seg = 1:2
        idx_cell{1,seg} = idx_imm_nonRew(ismember(idx_imm_nonRew, segment(:,seg))); 
        idx_cell{2,seg} = idx_mov_nonRew(ismember(idx_mov_nonRew, segment(:,seg)));
    end
    
    %% Divide rewards into early and late
    if ~isempty(idx_rew)
        segment = 3; % How many segments to divide recording into
        segment = floor(length(beh(x).reward)/segment); % How many rewards included in each segment
        rew_seg = beh(x).reward(1:segment); % EARLY reward delivery
        rewYes = extractRewardedTrials(rew_seg./beh(x).Fs, beh(x).lick./beh(x).Fs, [0 0.5]); % Extract rewarded trials
        rew_seg = rew_seg(rewYes); % EARLY
        idx_cell{3,1} = extractEventST([1:length(fp_mat)]', floor(rew_seg), floor(rew_seg)+50, 1);
        rew_seg = beh(x).reward(length(beh(x).reward)-segment+1:end); % LATE reward delivery
        rewYes = extractRewardedTrials(rew_seg./beh(x).Fs, beh(x).lick./beh(x).Fs, [0 0.5]); % Extract rewarded trials
        rew_seg = rew_seg(rewYes); % LATE
        idx_cell{3,2} = extractEventST([1:length(fp_mat)]', floor(rew_seg), floor(rew_seg)+50, 1);
    end

    %% Cross correlation, early or late motivational state
    mat(x).rec = beh(x).rec;
    fp_sub = []; 
    for z = 1:3
        for seg = 1:2 % EARLY or LATE motivational state
            if length(idx_cell{z,seg})< 2; continue; end
            fp_sub = fp_mat(idx_cell{z,seg},:); % signal during EARLY motivational state
            [corr_tmp, lags] = xcorr(fp_sub(:,1), fp_sub(:,2), 5*Fs, 'coeff'); % cross-correlation
            mat(x).corr{z} = corr_tmp;
            corr_cell{z,seg}(:,x) = corr_tmp;       % cross-correlation
        end
        fp_sub = [fp_mat(idx_cell{z,1},:); fp_mat(idx_cell{z,2},:)];
        fp_sub_new = fp_sub(:,2);
        tmp_shuff = []; 
        for s = 1:50
            fp_sub_new = circshift(fp_sub_new, Fs);
            tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_sub_new, 5*Fs, 'coeff');
        end
        corr_shuff{1,z}(:,x) = prctile(tmp_shuff, 5, 2); % shuffle 5th percentile
        corr_shuff{2,z}(:,x) = prctile(tmp_shuff, 50, 2); % shuffle 50th percentile
        corr_shuff{3,z}(:,x) = prctile(tmp_shuff, 95, 2); % shuffle 95th percentile
    end
    %% Shuffled cross-correlation
%     fp_sub_new = fp_mat(:,2);
%     tmp_shuff = []; 
%     for s = 1:50
%         fp_sub_new = circshift(fp_sub_new, Fs);
%         tmp_shuff(:,s) = xcorr(fp_mat(:,1), fp_sub_new, 5*Fs, 'coeff');
%     end
%     mat(x).shuff = prctile(tmp_shuff, [5 50 95], 2);
%     corr_shuff{1,1}(:,x) = prctile(tmp_shuff, 5, 2); % shuffle 5th percentile
%     corr_shuff{2,1}(:,x) = prctile(tmp_shuff, 50, 2); % shuffle 50th percentile
%     corr_shuff{3,1}(:,x) = prctile(tmp_shuff, 95, 2); % shuffle 95th percentile
        
%%
    waitbar(x/length(beh),h);
end
close(h);

%% N = X mice
corr_an = cell(3,2); corr_shuff_an = cell(3,1);
min_val = cell(1,2); min_lag = cell(1,2);

tmp = {}; for x = 1:length(beh); tmp{x} = strtok(beh(x).rec,'-'); end
uni = unique(tmp); nAn = length(uni);

for x = 1:nAn
    idx = strcmp(tmp,uni{x});
    for a = 1:3
        for b = 1:2
            corr_adj = corr_cell{a,b};
            corr_adj = corr_adj - nanmean(corr_adj([find(lags./Fs == -3):find(lags./Fs == -1)],:));
            corr_an{a,b}(:,x) = nanmean(corr_adj(:,idx),2);
        end
        corr_adj = corr_shuff{a,1};
        corr_shuff_an{a,1}(:,x) = nanmean(corr_adj(:,idx),2);
    end
end
for z = 1:2 % Iterate over motivational state
    for y = 1:3 % Iterate over behavioral states
        [min_val{z}(:,y), ii] = min(corr_an{y,z});
        min_lag{z}(:,y) = lags(ii)./50;
    end
    min_lag{z}(isnan(min_val{z})) = nan;
end