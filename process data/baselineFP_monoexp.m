function [dF, baseline] = baselineFP_exp(FP, Fs, baseWin)
% Baseline adjust photometry signal using mono or double exponential decay 
% function fit to a specified baseline portion of signal
% 
% [dF, baseline] = baselineFP_exp(FP, Fs, baseWin)
%
% INPUT:
% 'FP' - photometry signal to baseline
% 'Fs' - sampling frequency
% 'baseWin' - set baseline range from 1:X (in seconds)
%
% OUTPUT:
% 'dF_F' - baseline adjusted trace (%dF_F)
% 'baseline' - optional ouput of trend used to baseline
%
% Anne Krok, January 2026
%

% generate time vector for processing
time    = linspace(0, length(FP)/Fs, length(FP))';

% parse baseline portion of signal 
[~,idx_base] = min(abs(time - baseWin)); % baseline is first X seconds of trace
idx_base = 1:idx_base;
t_base  = time(idx_base); % baseline time vector
y_base  = FP(idx_base);   % baseline signal

% model: mono-exponential decay function, a*exp(-x/tau) + c
% ft = fittype('a*exp(-x/tau)+c','independent','x','coefficients',{'a','tau','c'});
% opts = fitoptions(ft);
% opts.StartPoint = [max(FP)-median(FP), (time(end)-time(1))/5, median(FP)];
% opts.Lower = [0, 0, -Inf];
% [fitted, ~] = fit(t_base, y_base, ft, opts);
% baseline = fitted.a*exp(-time./fitted.tau) + fitted.c;

% alternative model: double exponential decay function, a*exp(-x/tau1) + b*exp(-x/tau2) + c
ft = fittype('exp2');
opts = fitoptions(ft);
[fitted, ~] = fit(t_base, y_base, ft, opts);
baseline = feval(fitted, time);

% baseline-corrected signal (dF/F)
dF = (FP - baseline)./baseline;
dF = dF*100; % convert to %

end

%% Plot check
% figure;
% subplot(4,1,1); plot(time,FP); title('Raw signal'); ylabel('y');
% subplot(4,1,2); plot(time,baseline,'r-'); hold on; plot(t_base,y_base,'ko'); title('Fitted baseline trend'); ylabel('trend');
% subplot(4,1,3); plot(time,FP-baseline); title('De-bleached signal'); xlabel('time'); ylabel('y_{debleached}');
% subplot(4,1,4); plot(time,dF); title('dF/F'); xlabel('time'); ylabel('FP (dF/F)');