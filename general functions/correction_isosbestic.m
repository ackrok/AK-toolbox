function correctedsignal=correction_isosbestic(signal,isosbestic,fs)
%% create model from isosbestic channel

%getting rid of artifacts
isosbestic(1:round(fs*60),2)=isosbestic(round(fs*60),2);
signal(1:round(fs*60),2)=signal(round(fs*60),2);

isosbesticmodel=fit(isosbestic(:,1),isosbestic(:,2),'exp2');
%% trim to make sure they are the same length
if length(signal)>length(isosbestic)
    signal=signal(1:length(isosbestic),:);
elseif length(isosbestic)>length(signal)
        isosbestic=isosbestic(1:length(signal),:);
end

%% fit to signal
newmodel=isosbesticmodel(signal(:,1));

if newmodel(1) > signal(1,2)
    scale=newmodel(1)-signal(1,2);
    newY=newmodel-scale;
elseif signal(1,2)>newmodel(1)
    scale=signal(1,2)-newmodel(1);  
    newY=newmodel+scale;
end

%% getting corrected signal
newsignal=signal(:,2)-newY;
correctedsignal=interpbaseline(newsignal,'linear',round(fs*60));

%% dFF
baselineidx=find(signal(:,1)<=10*60);
dFF=(newsignal-mean(newsignal(baselineidx)))/mean((newsignal(baselineidx)));
zscoreddFF=(dFF-mean(dFF(baselineidx)))/std(dFF(baselineidx));

figure
plot(correctedsignal)
hold on
plot(newY)
plot(signal(:,2))
plot(isosbestic(:,2))
plot(isosbesticmodel(signal(:,1)))
title('new signal')
legend('corrected signal','scaled isosbestic model','original','isosbestic','isosbestic model')

figure
subplot(2,1,1)
plot(dFF)
title('dFF')
subplot(2,1,2)
plot(zscoreddFF)
title('z-score dFF')
