function [corrected_signal,dFF]=correction_red(signal,fs)

t     = signal(:,1);
F_raw = signal(:,2);

F_raw(1:60)=F_raw(60);


drugfit=fit(t,F_raw,'exp2'); %exponential fit to rGRABDA MDMA data
fittedmodel=feval(drugfit,t); %fits model onto data. This is predicted y based on time x
corrected_signal=F_raw-fittedmodel;
% corrected_signal=interpbaseline(newsignal,'linear',1800);

baseline=1:10*60*fs;
F0 = mean(corrected_signal(baseline)); 
% If F0 is too close to zero or negative, add a constant offset before dF/F:
if F0 <= 0 || abs(F0) < 1e-6
    offset = mean(F_raw(baseline));             % use raw baseline level as DC offset
    corrected_shifted = corrected_signal + offset;
    F0 = mean(corrected_shifted(baseline));
    dFF = (corrected_shifted - F0) / F0;
else
    dFF = (corrected_signal - F0) / F0;
end

% figure
% plot(corrected_signal)
% hold on
% plot(F_raw)
% hold on
% plot(ypdrug)
% legend('corrected signal','original signal','model')

figure;
subplot(3,1,1);
plot(t, F_raw, 'k'); hold on;
plot(t, fittedmodel, 'r');
xlabel('Time (s)'); ylabel('F');
legend('Raw', 'Stretched smoothed baseline');
title('Raw signal and stretched smoothed baseline');

subplot(3,1,2);
plot(t, corrected_signal, 'b');
xlabel('Time (s)'); ylabel('Corrected F');
title('Corrected signal (raw - stretched baseline)');

subplot(3,1,3);
plot(t, dFF, 'm');
xlabel('Time (s)'); ylabel('\DeltaF/F');
title('\DeltaF/F (using corrected signal)');

end