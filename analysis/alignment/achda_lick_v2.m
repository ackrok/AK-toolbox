% good_rew = [10:16,20,22:26,28:29,32:35,40:42,44:46];
% beh = modAChDA(good_rew);

%% Extract data and process lick and reward trials
out = struct;
window = [0 1]; % Window for reward to be collected

for x = 1:length(beh)
    
    out(x).rec = beh(x).rec;
    out(x).FP = beh(x).FP;
    out(x).FPnames = beh(x).FPnames;

%% Extract licks, correct for timing
    lick = beh(x).lick(:)/beh(x).Fs; % Licks, in seconds
    lick_repeat = [diff(lick.*1000) > 50]; % Identify licks that are <50ms after previous lick
    lick_sub = lick; lick_sub(1) = [];
    lick = [lick(1); lick_sub(lick_repeat)]; % Overwrite lick vector
    out(x).lick_new = lick; % Corrected licks, in seconds
    
%% Identify rewarded trials
rew = beh(x).reward(:)/beh(x).Fs; % Extract reward delivery times, adjusting event times to be in seconds
    out(x).delivery = rew; % Reward delivery times, in seconds

    rew_lick0 = []; rew_prewlick = [];
    bin = 1/1000; 
    peth = getClusterPETH(lick, rew, bin, window); % PETH: lick aligned to reward in 1 ms bins
    cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
    rew_lick0 = find(sum(cts,1) == 0); % Find reward index where total licks within window is 0

    % bin = 1/1000; window = [-0.25 0];
    % peth = getClusterPETH(lick, rew, bin, window);
    % rew_prewlick = find(sum(peth.cts{1},1) >= 1); % Find reward index for trials where mouse licks preceding reward

    rew([rew_lick0, rew_prewlick]) = nan; 
    cts(:, [rew_lick0, rew_prewlick]) = nan; % Remove non-rewarded trials and trials where mouse licks preceding reward
    ev = rew(~isnan(rew)); % Event time is rewarded trials
    
    out(x).rew_no = find(isnan(rew)); % Index of deliveries where animal did not lick to receive reward
    out(x).rew_yes = find(~isnan(rew)); % Index of delivieries where animal DID lick

%% Onset of lick bouts for each rewarded trial
    bin = 1/1000;
    peth = getClusterPETH(lick, ev, bin, window); % PETH: lick aligned to rewarded trials in 1 ms bins
    cts = peth.cts{1}; % Lick counts in 1ms bins for each reward trial
    [~, lick_first] = max(cts~=0, [], 1); % Find first non-zero index for each trial
%     [~, lick_last] = max(flipud(cts)~=0, [], 1); % Find last non-zero index
%     lick_last = window(2)/bin - (lick_last - 1); % Flip back because last lick was determined from flipped matrix
    out(x).rew_lick = lick_first(:).*bin + ev(:); % First lick after delivery for rewarded trials
    tmp = nan(length(out(x).delivery),1); % Initialize output vector
    tmp(out(x).rew_yes) = out(x).rew_lick;
    out(x).rew_lick = tmp; % Include NaNs for non-rewarded trials
    
%% Number of licks per reward period (0:2s)
    peth = getClusterPETH(lick, ev, 2, [0 2]); % PETH: lick aligned to rewarded trials
    cts = peth.cts{1}(:); % Lick counts for each reward trial across 2 second period post-delivery
    out(x).lick_num = nan(length(out(x).delivery),1);
    out(x).lick_num(out(x).rew_yes) = cts;
    
end
fprintf('Done re-processing lick and reward trials. Window: %d to %1.2f seconds \n',window(1),window(2));

%% Align DA and ACh to (a) rewarded trials, (b) non-rewarded trials, (c) onset of lick bout
align_rewYes = cell(length(out),2);
align_rewNo = cell(length(out),2);
align_lickOn = cell(length(out),2);
win = [-1 2]; Fs = 50;
for x = 1:length(out) % For each recording...
    for y = 1:length(out(x).FP)
        fp = out(x).FP{y}; % Extract photometry signal
        [sta, time] = getSTA(fp, out(x).delivery(out(x).rew_yes), Fs, win); % Align DA and ACh to delivery for rewarded trials
        align_rewYes{x,y} = sta; % Store into output cell array
        [sta] = getSTA(fp, out(x).delivery(out(x).rew_no), Fs, win); % Align DA and ACh to delivery for non-rewarded trials
        align_rewNo{x,y} = sta; % Store into output cell array
        [sta] = getSTA(fp, out(x).rew_lick, Fs, win); % Align DA and ACh to onset of lick bout
        align_lickOn{x,y} = sta; % Store into output cell array
    end
end

%% For each animal
tmp = {}; for x = 1:length(out); tmp{x} = strtok(out(x).rec,'-'); end 
uni = unique(tmp); % Unique animal ID
lbl = {'rewYes','rewNo','lickOn'};
align_uniAvg = cell(length(lbl),2); align_uniN = align_uniAvg;
for x = 1:length(uni) % For each animal...
    for y = 1:2
        ii = find(strcmp(tmp, uni{x})); % Identify which recordings correspond to this animal
        pull = horzcat(align_rewYes{ii,y}); % Concatenate across all recordings from one animal
        align_uniAvg{1,y}(:,x) = nanmean(pull,2); % Average across all
        align_uniN{1,y}(x) = size(pull,2);
        pull = horzcat(align_rewNo{ii,y}); 
        align_uniAvg{2,y}(:,x) = nanmean(pull,2);
        align_uniN{2,y}(x) = size(pull,2);
        pull = horzcat(align_lickOn{ii,y}); 
        align_uniAvg{3,y}(:,x) = nanmean(pull,2);
        align_uniN{3,y}(x) = size(pull,2);
    end
end
