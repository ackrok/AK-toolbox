function [corrected_signal,dFF]=correction_baseline(signal,fs)
% 
% %% Get 10 minute baseline
% baseline=find(signal(:,1)>0 & signal(:,1)<=10*60);
% signal(1:round(fs*60),2)=signal(round(fs*60),2);
% 
% 
% %% model baseline 
% 
% fittedmodel=resample(signal(baseline,2),length(signal(:,1)),baseline(end))';
% newsignal=signal(:,2)-fittedmodel;
% interpfactor=round(fs*60);
% corrected_signal=interpbaseline(newsignal,'linear',interpfactor);
% 
% 
% figure
% plot(corrected_signal)
% hold on
% plot(fittedmodel)
% hold on
% plot(signal(:,2))
% legend('corrected signal','fitted model','signal')
 %% new code
% signal: [N x 2], col1 = time (s), col2 = fluorescence
% fs    : sampling rate (Hz)

t     = signal(:,1);
F_raw = signal(:,2);

%% 1) Clamp first minute (artifact removal)
artifact_sec      = 2;                         % 1 minute
artifact_end_idx  = min(round(artifact_sec*fs), numel(F_raw));
F_raw(1:artifact_end_idx) = F_raw(artifact_end_idx);


%% 2) Define baseline: first 10 minutes (pre‑drug)
baseline_sec = 10*60;                           % 600 seconds
baseline     = t >= 0 & t <= baseline_sec;
F_base       = F_raw(baseline);


%% 3) Smooth the baseline so we keep only slow bleaching dynamics
%    Choose a fairly long window, e.g. 30–60 s, depending on your data
smooth_win_sec    = 60;                         % adjust if needed
smooth_win_samples = max(1, round(smooth_win_sec * fs));
F_base_smooth     = movmean(F_base, smooth_win_samples);

%% 4) Stretch (resample) smoothed baseline to full length
N_total    = numel(F_raw);
N_baseline = numel(F_base_smooth);

% This creates a stretched version of the baseline trend
fittedmodel = resample(F_base_smooth, N_total, N_baseline);
fittedmodel = fittedmodel(:);                   % make sure column

% Optional: soften edges a little to avoid step at start/end
edge_win = round(60*fs);                         % e.g. 5 s of smoothing
fittedmodel = movmean(fittedmodel, max(1,edge_win));
fittedmodel = fittedmodel(:);                 % make sure column

% Note: fittedmodel is now your "artificial bleaching curve" with same length as F_raw.
% It has the same value distribution as the smoothed baseline, just stretched in time.

%% 5) Subtract the stretched bleaching model from the original signal
corrected_signal = F_raw - fittedmodel;

%% 6) Compute ΔF/F using corrected signal and baseline only (optional)
%    Because we subtracted a baseline-like curve, the corrected baseline
%    might be near zero; if ΔF/F is unstable, add a constant offset as shown below.

F0 = mean(corrected_signal(baseline));          % might be small; inspect it

% If F0 is too close to zero or negative, add a constant offset before dF/F:
if F0 <= 0 || abs(F0) < 1e-6
    offset = mean(F_raw(baseline));             % use raw baseline level as DC offset
    corrected_shifted = corrected_signal + offset;
    F0 = mean(corrected_shifted(baseline));
    dFF = (corrected_shifted - F0) / F0;
else
    dFF = (corrected_signal - F0) / F0;
end

%% 7) Quick diagnostics
fprintf('Baseline mean F_raw = %.4g, mean fittedmodel(baseline) = %.4g, F0 = %.4g\n', ...
        mean(F_raw(baseline)), mean(fittedmodel(baseline)), F0);

%% 8) Plot for visual inspection
figure;
subplot(3,1,1);
plot(t, F_raw, 'k'); hold on;
plot(t, fittedmodel, 'r');
xlabel('Time (s)'); ylabel('F');
xline(600,'--')
legend('Raw', 'Stretched smoothed baseline');
title('Raw signal and stretched smoothed baseline');

subplot(3,1,2);
plot(t, corrected_signal, 'b');
xlabel('Time (s)'); ylabel('Corrected F');
title('Corrected signal (raw - stretched baseline)');
xline(600,'--')

subplot(3,1,3);
plot(t, dFF, 'm');
xlabel('Time (s)'); ylabel('\DeltaF/F');
title('\DeltaF/F (using corrected signal)');
xline(600,'--')

end