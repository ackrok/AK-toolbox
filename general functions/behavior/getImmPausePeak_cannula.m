function [cannula] = getImmPausePeak_cannula(cannula)
%Extract characteristics of troughs and peaks of ACh signal during periods
%of immobility, for cannula infusion experiments
%
%   [cannula] = getImmPausePeak_cannula(cannula)
%
%   Description: This function  analyzing fluctuations of ACh signal
%   during periods of immobility, and extracting the frequency, duration,
%   and amplitude of the troughs and peaks of fluctuations.
%       Frequency is defined as the number of events per second.
%       Duration is defined as the width at the half maximum amplitude
%       Amplitude is defined as the amplitude from zero
%
%   INPUT
%   'cannula' - structure containing photometry, behavior, infusion
%       information for multiple recordings, using extractCannula
%
%   OUTPUT
%   'cannula' - updated structure
%
% Anya Krok, September 2022

%% INPUTS
NumStd = 1.5; % Number of standard deviations of filtered photometry
% signal during immobility. Troughs and peaks are identified that
% have an amplitude that is more than specified value * STD
s = cannula; % renaming structure to shorter name

%%
params = [];
for y = 1:length(s)
    %%
    beh = s(y).s; % extract beh sub-structure within cannula structure
    amp = []; dur = []; freq = []; % initialize output matrix for this infusion
%     ach2ach = cell(length(beh),2); % initialize output cell array for this infusion
%     da2ach = cell(length(beh),2);

    for x = 1:length(beh) % iterate over animal
        if isempty(beh(x).Fs); continue; end
        Fs = beh(x).Fs; % sampling frequency
        fp_mat = beh(x).FP{1};
        fp_mat = fp_mat - nanmean(fp_mat);
        fp_mat(:,2) = beh(x).FP{2} - nanmean(beh(x).FP{2}); 
        
        idx_inf = [s(y).win(x,1).*(Fs*60) : s(y).win(x,2).*(Fs*60)]'; % infusion window
        rewWindow = 2*Fs; % how many samples after reward delivery is the reward window
        idx_inf_rew = extractEventST(idx_inf, floor(beh(x).reward), floor(beh(x).reward)+rewWindow, 1); % identify sample during reward
        idx_inf_rest = extractEventST(idx_inf, beh(x).onRest, beh(x).offRest, 1); % index: infusion + rest
        idx_inf_rest = idx_inf_rest(~ismember(idx_inf_rest, idx_inf_rew)); % exclude reward, include rest

        %% Bandpass filter
        Fpass = [0.1 10];
        Fs = 50; %sampling rate, has to be at least double of your high pass frequency
        Wn = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; 
        [b,a] = butter(3,Wn);
        data_filt = filtfilt(b,a,fp_mat(:,1)); % signal is your photometry data, output is the filtered data

        %% Instantaneous phase
        data_filt = data_filt(idx_inf_rest); % IMMOBILITY ONLY
        H = hilbert(double(data_filt));
        data_phase  = angle(H); % output is the instantaneous phase
        fp_phase = data_phase;
        fp_deg = rad2deg(data_phase);

        %Take the absolute of the filtered signal and calculate the standard deviation
        rmssig  = abs(data_filt);
        stdsig = std(rmssig);        

        %Find the index of peaks that are more than NumStd standard deviations (NumStd could be 1.5 standard deviations)
        % finds the 0 degree phase indices that cross 1.5 standard deviations. Might have to remove the ; depending on whether your data is a row or column array
        if y == 1; params(x) = NumStd*stdsig; end
        idxPeak = find(data_filt>(params(x)) & [0; diff(data_phase>0)]); 
        idxPause = find(-data_filt>(params(x)) & [0; diff(-data_phase>0)]); 
%         figure; hold on
%         plot(data_filt, 'k')
%         stem(idxPeak, 10*ones(length(idxPeak),1), 'g'); 
%         stem(idxPause, 10*ones(length(idxPause),1), 'r');

        %% Peak characterization
        tmp_amp_peak = data_filt(idxPeak); 
        tmp_amp_pause = data_filt(idxPause);
        tmp_dur_peak = []; tmp_dur_pause = []; win = 1*Fs;
        for z = 1:length(idxPeak)
            halfMax = 0.5*data_filt(idxPeak(z)); % amplitude at half-max
            if idxPeak(z)-win < 0
                a = [data_filt(1 : idxPeak(z))]; a = flipud(a); % if peak is close to start of recording
            else
                a = [data_filt(idxPeak(z)-win : idxPeak(z))]; a = flipud(a); % segment preceding idx of maximum deflection
            end
            b = find(a < 0, 1, 'first') - 1; 
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            c = []; [~,c(1)] = min(abs(a-halfMax)); c(1) = c(1)-1; % idx at half max
            if idxPeak(z)+win > length(data_filt)
                a = [data_filt(idxPeak(z) : end)]; 
            else
                a = [data_filt(idxPeak(z) : idxPeak(z)+win)]; % segment following idx of maximum deflection
            end
                
            b = find(a < 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            [~,c(2)] = min(abs(a-halfMax)); % idx at half max
            
            tmp_dur_peak(z) = sum(c);
        end % width at half max for peaks
        for z = 1:length(idxPause)
            halfMax = 0.5*data_filt(idxPause(z)); % amplitude at half-max
            if idxPause(z)-win < 0
                a = [data_filt(1 : idxPause(z))]; a = flipud(a); % if peak is close to start of recording
            else
                a = [data_filt(idxPause(z)-win : idxPause(z))]; a = flipud(a); % segment preceding idx of maximum deflection
            end
            b = find(a > 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            c = []; [~,c(1)] = min(abs(a-halfMax)); c(1) = c(1)-1; % idx at half max
            if idxPause(z)+win > length(data_filt)
                a = [data_filt(idxPause(z) : end)]; 
            else
                a = [data_filt(idxPause(z) : idxPause(z)+win)]; % segment following idx of maximum deflection
            end
            b = find(a > 0, 1, 'first') - 1;
            if ~isempty(b); a = a(1:b); end % retain only values above 0
            [~,c(2)] = min(abs(a-halfMax)); % idx at half max
            
            tmp_dur_pause(z) = sum(c);
        end % width at half max for pauses

        amp(x,1) = nanmean(tmp_amp_pause); % Amplitude of maximum deflection
        amp(x,2) = nanmean(tmp_amp_peak); 
        dur(x,1) = nanmean(tmp_dur_pause); % Duration at half max
        dur(x,2) = nanmean(tmp_dur_peak); 
        freq(x,1) = (1./nanmean(diff(idxPause)))*Fs; % Frequency of maximum deflections
        freq(x,2) = (1./nanmean(diff(idxPeak)))*Fs; % Frequency of maximum deflections
        
        %% Photometry to peaks/pauses
%         [sta, t] = getSTA(fp_mat(idx_inf_rest,1), idxPause/Fs, Fs, [-6 2]);
%         ach2ach{x,1} = sta;
%         sta = getSTA(fp_mat(idx_inf_rest,1), idxPeak/Fs, Fs, [-6 2]);
%         ach2ach{x,2} = sta;
%         sta = getSTA(fp_mat(idx_inf_rest,2), idxPause/Fs, Fs, [-6 2]);
%         da2ach{x,1} = sta;
%         sta = getSTA(fp_mat(idx_inf_rest,2), idxPeak/Fs, Fs, [-6 2]);
%         da2ach{x,2} = sta;
    end
    s(y).params = params;
    s(y).amp = amp; 
    s(y).dur = dur; 
    s(y).freq = freq;
%     s(y).ach2ach = ach2ach;
%     s(y).da2ach = da2ach;
end
cannula = s;
end
