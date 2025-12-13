%%
good_rew = [10:16,20,22:26,28:29,32:35,40:42,44:46];
beh = modAChDA(good_rew);

%%
mat = struct;
h = waitbar(0,'coherogram');
rewWin = 1; % seconds in reward window
Fs = 50; % sampling rate
for x = 1:length(beh); aa = [2 1];
  
    mat(x).rec = beh(x).rec; % load recording name
    fp_mat = [];
    fp_mat(:,1) = [beh(x).FP{aa(1)} - nanmean(beh(x).FP{aa(1)})];
    fp_mat(:,2) = [beh(x).FP{aa(2)} - nanmean(beh(x).FP{aa(2)})]; 
    
    %% Behavioral states
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+(rewWin*Fs), 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    
    %% Divide data to assess early/late motivational state
    c = cell(3,2);
    segment = 3; % how many segments to divide recording into
    segment = reshape([1:length(fp_mat)],[length(fp_mat)/segment, segment]);
    segment(:,2) = []; % remove center segment, keep only early and late
    for seg = 1:2
        c{1,seg} = idx_imm_nonRew(ismember(idx_imm_nonRew, segment(:,seg))); 
        c{2,seg} = idx_mov_nonRew(ismember(idx_mov_nonRew, segment(:,seg)));
    end
    
    %% Divide rewards into early and late
    if ~isempty(idx_rew)
        segment = 3; % How many segments to divide recording into
        segment = floor(length(beh(x).reward)/segment); % How many rewards included in each segment
        rew_seg = beh(x).reward(1:segment); % EARLY reward delivery
        rewYes = extractRewardedTrials(rew_seg./beh(x).Fs, beh(x).lick./beh(x).Fs, [0 0.5]); % Extract rewarded trials
        rew_seg = rew_seg(rewYes); % EARLY
        c{3,1} = extractEventST([1:length(fp_mat)]', floor(rew_seg), floor(rew_seg)+(rewWin*Fs), 1);
        rew_seg = beh(x).reward(length(beh(x).reward)-segment+1:end); % LATE reward delivery
        rewYes = extractRewardedTrials(rew_seg./beh(x).Fs, beh(x).lick./beh(x).Fs, [0 0.5]); % Extract rewarded trials
        rew_seg = rew_seg(rewYes); % LATE
        c{3,2} = extractEventST([1:length(fp_mat)]', floor(rew_seg), floor(rew_seg)+(rewWin*Fs), 1);
    end
    
    %%
    for y = 1:size(c,1) % iterate over behavioral states
        for z = 1:size(c,2) % iterate over segments
            if ~isempty(c{y,z})
                sig1 = fp_mat(c{y,z},1); % extract indexes samples 
                sig2 = fp_mat(c{y,z},2); % extract indexes samples 
                [coher,ph,t,f] = bz_MTCoherogram(sig1,sig2,'frequency',50,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
                mat(x).coher{z}(:,y) = nanmean(coher,2); % collapse time dimension
                mat(x).phase{z}(:,y) = nanmean(ph,2); % collapse time dimension
            else
                mat(x).coher{z}(:,y) = nan(103,1);
                mat(x).phase{z}(:,y) = nan(103,1);
            end
        end
    end
    %%
    tmp_coher = []; tmp_phase = []; 
    new_2 = fp_mat(:,2);
    for s = 1:50 % repeat shuffle N times
        % new_2 = circshift(new_2, Fs);
        new_2 = new_2(randperm(length(new_2)));
        [c,p] = bz_MTCoherogram(fp_mat(:,1),new_2,'frequency',50,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
        tmp_coher(:,s) = nanmean(c,2); % collapse time dimension
        tmp_phase(:,s) = nanmean(p,2); % collapse time dimension
    end
    mat(x).coher_shuff = prctile(tmp_coher, [5 50 95], 2); % 5th, 50th, 95th percentiles
    mat(x).phase_shuff = prctile(tmp_phase, [5 50 95], 2); % 5th, 50th, 95th percentiles
    
    waitbar(x/length(beh),h);
end
close(h); fprintf('Coherogram Analysis Done !\n');

%% Extract IMM MOV REW
coher_state = cell(3,2); phase_state = cell(3,2);
for x = 1:length(mat)
    for y = 1:3; for z = 1:2
        if y == 3 && size(mat(x).coher{z},2) < 3
            coher_state{y,z}(:,x) = nan(103,1); 
            phase_state{y,z}(:,x) = nan(103,1);
        else
        coher_state{y,z}(:,x) = mat(x).coher{z}(:,y);
        phase_state{y,z}(:,x) = mat(x).phase{z}(:,y);
        end
    end; end
end
fprintf('Coherogram Analysis: extracted behavioral states\n');

%%
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni);
coher_an = cell(3,2); phase_an = cell(3,2);
coher_shuff = cell(3,1); phase_shuff = cell(3,1);
for x = 1:nAn
    ii = find(strcmp(rec,uni{x})); % index of matching recordings for this animal
    coh_tmp = [mat(ii).coher_shuff];
    ph_tmp = [mat(ii).phase_shuff];
    for y = 1:3; for z = 1:2 % iterate over behavioral states; iterate over segments
        coher_an{y,z}(:,x) = nanmean(coher_state{y,z}(:,ii),2); % average recordings for each animal, for each behavioral state
        phase_an{y,z}(:,x) = nanmean(phase_state{y,z}(:,ii),2); % average recordings for each animal, for each behavioral state
        coher_shuff{y}(:,x) = nanmean(coh_tmp(:,[y:3:size(coh_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
        phase_shuff{y}(:,x) = nanmean(ph_tmp(:,[y:3:size(ph_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
        end; end
end
fprintf('Coherogram Analysis: averaged for each mouse\n');