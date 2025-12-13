function [WF] = getWFprop(WF, varargin)
%Charactirization of waveform shape and temporal properties
% 
% Description: This function will calculate waveform features for units in
%   a recording, using waveforms extracted into T-Lab data structure
%
% [WFprop] = getWFprop(WF, varargin)
%
% INPUT
%   WF - structure generated with function loadWF, includes:
%       'WF.topMu' - average waveform for each unit
%       'WF.time' - vector of time values for plotting waveforms
%   varargin
%       {1} 'units' - vector of units to iterate waveform analysis over
%
% OUTPUT
%   WF - structure updated to now include:
%       'WF.peakLatency' = (*ms*) time between WF minimum and second maximum. 
%       'WF.halfWidth' = (*ms*) time width of the spike at half-maximal amplitude of the first maximum. 
%       'WF.amplitude' = voltage difference between the second maximum and the minimum. 
%       'WF.asymmetry' = ratio of the amplitude of the second maximum to that amplitude of the first maximum.
%       'WF.paramIdx' = idx used to calculate WF properties
%       'WF.paramVal' = values at idx used to calculate WF properties
%
% Created By: Anya Krok, June 2019
% Updated On: January 2020 
%

%%
if nargin == 2
    units = varargin{1};
    if iscolumn(units); units = units'; end %make row vector
elseif nargin == 1
    units = 1:length(WF);
end

%% INITIALIZE VARIABLES
% Peak latency = time between the minimum and the second maximum. 
% Spike half-width = width of the spike at half-maximal amplitude of the first maximum. 
% Voltage amplitude = voltage difference between the second maximum and the minimum. 
% Asymmetry = ratio of the amplitude of the second maximum to that amplitude of the first maximum.
peakLatency = zeros(length(WF),1); 
halfWidth   = zeros(length(WF),1);
amplitude   = zeros(length(WF),1); 
asymmetry   = zeros(length(WF),1); 
paramIdx    = zeros(length(WF),5); 
paramVal    = zeros(length(WF),4);

%% AUTOMATIC PROCESSING
fprintf('waveform automatic processing - '); 
fig = figure; set(fig,'Position',get(0,'Screensize'))
plotM = floor(sqrt(length(units))); plotN = ceil(length(units)/plotM);

for x = 1:length(units) %CHANGE
    n = units(x);   %CHANGE
    waveform = WF(n).topMu; time = WF(n).time;
    %% Plot waveform for each unit
    sp(x) = subplot(plotM,plotN,x); hold on
        h1 = plot(waveform,'k'); %plot for ginput
        title(['unit #',num2str(n),' - waveform']); 
        ylabel('amplitude (uV)'); xlabel('index');
        
    %% Automatically select boundaries    
    [minY,minI] = min(waveform); %find index of minimum point    
    X = [1, minI-1, minI+1, length(waveform)]; %auto-boundaries    
    Y = waveform(X);    
    %h2 = scatter(X,Y,'filled','k');    
    
    %% Waveform properties
    [out] = waveformPropHelper(X, waveform, time);
        halfIdx = [out(5).val(2),out(5).val(4)]; 
        minIdx = out(5).val(3); minVal = out(6).val(3);
        maxIdx = [out(5).val(1),out(5).val(5)]; maxVal = [out(6).val(1),out(6).val(4)];
    h3 = scatter(halfIdx,waveform(halfIdx),'filled','r'); 
    h4 = scatter(minIdx,minVal,'filled','g');
    h5 = scatter(maxIdx,maxVal,'filled','b');
    xlim([X(1) X(4)])
    
    %% Load parameters for this unit        
    halfWidth   (n)   = out(1).val;
    peakLatency (n)   = out(2).val;
    amplitude   (n)   = out(3).val;
    asymmetry   (n)   = out(4).val;
    paramIdx    (n,:) = out(5).val;
    paramVal    (n,:) = out(6).val;
end
linkaxes(sp,'x')
fprintf('DONE! \n'); 

%% MANUAL ADJUSTMENT
answer = inputdlg({'Enter space-separated UNIT# for manual adjustment:','for deletion:'},...
    'Input',[1 60]);
adj = str2num(answer{1}); 
del = str2num(answer{2});

% clear previous generated parameters
halfWidth([adj,del]) = 0; peakLatency([adj,del]) = 0; 
amplitude([adj,del]) = 0; asymmetry([adj,del]) = 0;
paramIdx([adj,del],:) = 0; paramVal([adj,del],:) = 0;
close(fig);

%% MANUAL PROCESSING
if ~isempty(adj)
    fprintf('waveform manual processing - '); 
    h = waitbar(0,'processing waveforms');
    for x = 1:length(adj)
        n = adj(x);
        fprintf('unit %d. ', n);
        waveform = WF(n).topMu; time = WF(n).time;
        
        %% Manually select boundaries 
        fig = figure; hold on;
            set(fig,'Position',get(0,'Screensize'));
            h1 = plot(waveform,'k'); %plot for ginput
            title(['unit #',num2str(n),' - (4x) boundaries for 1st max, min, 2nd max']);
            ylabel('amplitude (uV)'); xlabel('index');
            
        choice = 1;
        while choice == 1 %continue adjusting window until peaks look correct
            [X,~] = ginput(4);  %manually select 4 boundaries for 1st max, min, 2nd max
            X(X>length(waveform)) = length(waveform); X(X<1) = 1;
            X = round(X); Y = waveform(X);    
            h2 = scatter(X,Y,'filled','k'); %plot selected points
            
            [out] = waveformPropHelper(X, waveform, time);
                halfIdx = [out(5).val(2),out(5).val(4)]; 
                minIdx = out(5).val(3); minVal = out(6).val(3);
                maxIdx = [out(5).val(1),out(5).val(5)]; maxVal = [out(6).val(1),out(6).val(4)];
            h3 = scatter(halfIdx,waveform(halfIdx),'filled','r'); %plot 
            h4 = scatter(minIdx,minVal,'filled','g');
            h5 = scatter(maxIdx,maxVal,'filled','b');
            xlim([1 length(waveform)])
            
            [choice] = menu(['WF boundaries good for unit(',num2str(n),')'],'GO','DONE','SKIP','EXIT');
            delete(h2); delete(h3); delete(h4); delete(h5);
        end
        
        switch choice
            case 3; close(fig); continue %skips rest of for loop and begins next iteration for next cluster
            case 4; close(fig); break %exit loop completely
        end
        pause(2); close all;
        
        %% Store values
        halfWidth   (n)   = out(1).val;
        peakLatency (n)   = out(2).val;
        amplitude   (n)   = out(3).val;
        asymmetry   (n)   = out(4).val;
        paramIdx    (n,:) = out(5).val;
        paramVal    (n,:) = out(6).val;
        %%
        waitbar(find(n == adj)/length(adj));
    end
    close(h); %close waitbar, close figure
    fprintf('DONE! \n');
end

%% SAVING WAVEFORM DATA TO STRUCTURE
temp = WF;
if nargin == 2 && isfield(temp,'paramIdx') 
    populate = units; %if already exist, then only merge currently processed units
else 
    populate = 1:length(temp); %if does not exist, then fill completely
end
for n = populate
    temp(n).halfWidth   = halfWidth(n); 
    temp(n).peakLatency = peakLatency(n); 
    temp(n).amplitude   = amplitude(n); 
    temp(n).asymmetry   = asymmetry(n); 
    temp(n).paramVal    = paramVal(n,:); 
    temp(n).paramIdx    = paramIdx(n,:); 
end
WF = temp;
end

function [out] = waveformPropHelper(x,waveform,time,varargin)
% [out] = waveformPropHelper(x,waveform,time)
% 
% If inverted waveform:
%   [out] = waveformPropHelper(x,waveform,time,1)
%
% Description: This function uses indices (X), time vector, and waveform
% voltage to calculate waveform properties (half-width, peak-latency,
% amplitude, asymmetry). Output will also return the indices and values
% used to calculate these parameters.
% To be run within getWFprop.m
% 
% Anya, January 2020
% 
%% Values from indices
    % Input vector 'X' is indices. We will use time and voltage to determine
    % the values (time, voltage) that correspond to indices in X.
        % 1st maximum - local maximum between 1st and 2nd index
        % Minimum - local minimum between 2nd and 3rd index
        % 2nd maximum - local maximum between 3rd and 4th index
    [maxWave1 maxIdx1]  = max(waveform([x(1):x(2)]));     maxIdx1 = x(1)+maxIdx1-1; %adjust for range
    [minWave minIdx]    = min(waveform([x(2):x(3)]));     minIdx  = x(2)+minIdx-1;
    [maxWave2, maxIdx2] = max(waveform([x(3):x(4)]));     maxIdx2 = x(3)+maxIdx2-1;
    
    % If waveform is inverted, then we will use 1st min , maximum, 2nd min
    if maxWave1 == minWave
        [maxWave1, maxIdx1] = min(waveform([x(1):x(2)]));  maxIdx1 = x(1)+maxIdx1-1;  
        [minWave, minIdx]   = max(waveform([x(2):x(3)]));  minIdx  = x(2)+minIdx-1;    
        [maxWave2, maxIdx2] = min(waveform([x(3):x(4)]));  maxIdx2 = x(3)+maxIdx2-1;   
    end
    
%% Determine half-width
    % Voltage at half-width is calculated based on 1st maximum and minimum
    halfMax = (minWave + maxWave1)/2;            
    % Then, find index along WFmean closest to half-width voltage
    % half-width idx should be between 1st max idx and 2nd max idx
    [~,halfIdx1] = min(abs(waveform(maxIdx1:minIdx)-halfMax));  halfIdx1 = halfIdx1+maxIdx1-1; %adjust for range
    [~,halfIdx2] = min(abs(waveform(minIdx:maxIdx2)-halfMax));  halfIdx2 = halfIdx2+minIdx-1;
    
%% Calculate parameters: half-width, peak-latency, amplitude, asymmetry
    % halfWidth (ms) - width of spike at half maximal amplitude
    % peakLatency (ms) - time from max deflection to second max
    % amplitude (uV) - voltage difference between max deflection and second max
    % asymmetry - ratio of amplitude of second max and first max
    % paramIdx = idx used 
    % paramVal = values at idx
    halfWidth   = time(halfIdx2) - time(halfIdx1);
    peakLatency = time(maxIdx2)  - time(minIdx); 
    amplitude   = maxWave2 - minWave;
    asymmetry   = maxWave2/maxWave1;
    paramIdx    = [maxIdx1, halfIdx1, minIdx, halfIdx2, maxIdx2];
    paramVal    = [maxWave1, halfMax, minWave, maxWave2];
    
%% Output variable
    out(1).name = 'halfWidth';      out(1).val = halfWidth;
    out(2).name = 'peakLatency';    out(2).val = peakLatency;
    out(3).name = 'amplitude';      out(3).val = amplitude;
    out(4).name = 'asymmetry';      out(4).val = asymmetry;
    out(5).name = 'paramIdx';       out(5).val = paramIdx;
    out(6).name = 'paramVal';       out(6).val = paramVal;

end
