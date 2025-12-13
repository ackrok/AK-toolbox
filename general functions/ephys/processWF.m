function [WF] = processWF(WF, varargin)
%Process spike waveforms
% 
% Description: This function will process spike waveforms for units and 
%   compute the coordinates of waveform features to be used to compute
%   waveform properties 
%
% [WF] = processWF(WF, varargin)
%
% INPUT
%   'WF' - structure with extracted waveforms generated using loadWF
%   varargin:
%       {1} 'units' - vector of units to iterate waveform analysis over
%
% OUTPUT
%   'WF' - structure updated to include coordinates for waveform features:
%       'WF.ptX' = x-values, in ms, for waveform features
%           ptX = [max_1, min, max_2, half_1, half_2]
%
%       'WF.ptY' = y-values, in uV, for waveform features
%           ptY = [max_1, min, max_2, halfmax]
%
% Created By: Anya Krok, June 2019
% Updated On: January 2020 - added automatic processing
% Updated On: March 2020 - split getWFprop into two functions: processWF and getWFprop
%

%%
switch nargin
    case 1
        units = [1:length(WF)]; %Will process all units
    case 2
        units = varargin{1}; %Will process only specified units
        if iscolumn(units)
            units = units'; %Make into row vector
        end
end

%% Initialize Variables
xVal = zeros(length(WF),5); 
yVal = zeros(length(WF),4);

%% AUTOMATIC PROCESSING
fprintf('\n Waveform automatic processing - '); 
fig = figure; 
fig.Position = [1449 754 2545 1303]; % set(fig,'Position',get(0,'Screensize'))
plotM = floor(sqrt(length(units))); plotN = ceil(length(units)/plotM);

for z = 1:length(units) 
    n = units(z);   
    waveform = WF(n).topMu; %Extract this unit's waveform from WF structure
    time = WF(n).time; 
    %% Plot waveform for each unit
    sp(z) = subplot(plotM,plotN,z); hold on
        h1 = plot(waveform,'k'); %Plot waveform
        title(['unit #',num2str(n),' - waveform']); 
        xlabel('index'); ylabel('amplitude (uV)'); 
        xlim([1 length(waveform)]); 
        
    %% Automatically select boundaries    
    [~,minI] = min(waveform); %Find index of point of maximum deflection
    if minI == 1 || minI == length(waveform)
        continue
    end
    X = [1, minI-1, minI+1, length(waveform)]; %Auto-set feature boundaries    
    Y = waveform(X);    
    %h2 = scatter(X,Y,'filled','k');    
    
    %% Waveform features
    [ptIdx, ptVal] = waveformPropHelper(X, waveform);
    % ptIdx = [max1idx, minIdx, max2idx, half1idx, half2idx]
    % ptVal = [max1val, minVal, max2val, halfMax]
    maxIdx  = [ptIdx(1), ptIdx(3)];     maxVal  = [ptVal(1), ptVal(3)];
    minIdx  = ptIdx(2);                 minVal  = ptVal(2); 
    halfIdx = [ptIdx(4), ptIdx(5)];     halfVal = waveform(halfIdx);
    %Plot computed waveform feature index and value overlying waveform    
    h3 = scatter(halfIdx, halfVal, 'filled', 'r'); 
    h4 = scatter(minIdx, minVal, 'filled', 'g');
    h5 = scatter(maxIdx, maxVal, 'filled', 'b');
    xlim([X(1) X(4)])
    
    %% Store values
    xVal (n,:) = time(ptIdx); %x-values, in ms, for waveform features
    yVal (n,:) = ptVal; %y-values, in uV, for waveform features
end
linkaxes(sp,'x')
fprintf('DONE! \n'); 

%% MANUAL ADJUSTMENT
answer = inputdlg({'Enter space-separated UNIT# for manual adjustment:','for deletion:'},...
    'Input',[1 60]);
adj = str2num(answer{1}); 
del = str2num(answer{2});

% clear values generated during automatic processing
xVal([adj,del],:) = 0; 
yVal([adj,del],:) = 0;
close(fig);

%% MANUAL PROCESSING
if ~isempty(adj)
    fprintf('Waveform manual processing - '); 
    h = waitbar(0,'Manually Processing Waveforms');
    for z = 1:length(adj)
        n = adj(z);
        fprintf('unit %d. ', n);
        waveform = WF(n).topMu; time = WF(n).time;
        
        %% Manually select boundaries 
        fig = figure; hold on;
            fig.Position = [1 41 2560 1323]; %set(fig,'Position',get(0,'Screensize'));
            h1 = plot(waveform,'k'); %Plot waveform, will be used for ginput
            title(['unit #',num2str(n),' - (4x) boundaries for 1st max, min, 2nd max']);
            ylabel('amplitude (uV)'); xlabel('index'); 
            xlim([1 length(waveform)]); 
            
        [choice] = menu(['WF boundaries good for unit(',num2str(n),')'],'GO','DONE','SKIP','EXIT');
        while choice == 1 %continue adjusting window until peaks look correct
            [X,~] = ginput(4);  %manually select 4 boundaries for 1st max, min, 2nd max
            X(X>length(waveform)) = length(waveform); X(X<1) = 1;
            X = round(X); Y = waveform(X);    
            h2 = scatter(X,Y,'filled','k'); %plot selected points
            
            [ptIdx, ptVal] = waveformPropHelper(X, waveform);
            % ptIdx = [max1idx, minIdx, max2idx, half1idx, half2idx]
            % ptVal = [max1val, minVal, max2val, halfMax]
            maxIdx  = [ptIdx(1), ptIdx(3)];     maxVal  = [ptVal(1), ptVal(3)];
            minIdx  = ptIdx(2);                 minVal  = ptVal(2); 
            halfIdx = [ptIdx(4), ptIdx(5)];     halfVal = waveform(halfIdx);
            %Plot computed waveform feature index and value overlying waveform    
            h3 = scatter(halfIdx, halfVal, 'filled', 'r'); 
            h4 = scatter(minIdx, minVal, 'filled', 'g');
            h5 = scatter(maxIdx, maxVal, 'filled', 'b');
            xlim([1 length(waveform)]); 
            
            [choice] = menu(['WF boundaries good for unit(',num2str(n),')'],'GO','DONE','SKIP','EXIT');
            delete(h2); delete(h3); delete(h4); delete(h5);
        end
        
        switch choice
            case 3; close(fig); continue %skips rest of for loop and begins next iteration for next cluster
            case 4; close(fig); 
                fprintf('\n Manual Processing NOT completed for units: ['); fprintf(' %d',adj(z:end)); fprintf(']\n');
                break %exit loop completely
        end
        pause(2); close all;
        
        %% Store values
        xVal (n,:) = time(ptIdx); %x-values, in ms, for waveform features
        yVal (n,:) = ptVal; %y-values, in uV, for waveform features
        %%
        waitbar(find(n == adj)/length(adj));
    end
    close(h); %close waitbar, close figure
    fprintf('DONE! \n');
end

%% SAVING WAVEFORM DATA TO STRUCTURE
temp = WF;
if nargin == 2 && isfield(temp,'ptX') 
    populate = units; %if already exist, then only over-write currently processed units
else 
    populate = 1:length(temp); %if does not exist, then fill completely
end
for n = populate
    temp(n).ptX = xVal(n,:); %x-values, in ms, for waveform features
    temp(n).ptY = yVal(n,:); %y-values, in uV, for waveform features
end
WF = temp;
end

function [ptIdx, ptVal] = waveformPropHelper(xIn, waveform)
% [ptIdx, ptVal] = waveformPropHelper(xIn, waveform)
%
% Description: This function uses indices, xIn, and unit waveform vector.
% Returns index and value for waveform features:
%   (a) first maximum
%   (b) point of maximum deflection
%   (c) second maximum
%   (d) spike half-maximums.
%
% OUTPUT
%   'ptIdx': ptIdx = [max_1_idx, min_idx, max_2_idx, half_1_idx, half_2_idx]
%   'ptVal': ptVal = [max_1_val, min_val, max_2_val, halfmax]
% 
% Created by: Anya Krok, January 2020
% 
%% Values from indices
    % Input vector 'X' is indices. We will use waveform vector to determine
    % the values (time, voltage) that correspond to indices in X.
        % 1st maximum - local maximum between 1st and 2nd index
        % Minimum   - local minimum between 2nd and 3rd index
        % 2nd maximum - local maximum between 3rd and 4th index
    [max_1_val, max_1_idx]  = max(waveform([xIn(1):xIn(2)]));   max_1_idx = xIn(1)+max_1_idx-1; %adjust for range
    [min_val, min_idx]    = min(waveform([xIn(2):xIn(3)]));     min_idx  = xIn(2)+min_idx-1;
    [max_2_val, max_2_idx]  = max(waveform([xIn(3):xIn(4)]));   max_2_idx = xIn(3)+max_2_idx-1;
    
    % If waveform is inverted, then we will use 1st min, maximum, 2nd min
    if max_1_val == min_val
        [max_1_val, max_1_idx] = min(waveform([xIn(1):xIn(2)]));  max_1_idx = xIn(1)+max_1_idx-1;  
        [min_val, min_idx]   = max(waveform([xIn(2):xIn(3)]));    min_idx  = xIn(2)+min_idx-1;    
        [max_2_val, max_2_idx] = min(waveform([xIn(3):xIn(4)]));  max_2_idx = xIn(3)+max_2_idx-1;   
    end
    
%% Determine half-width
    % Voltage at half-width is calculated based on 1st maximum and minimum
    halfmax = (min_val + max_1_val)/2;            
    % Then, find index by computing values closest to half-width
    % half-width idx should be between 1st max idx and 2nd max idx
    [~,half_1_idx] = min(abs(waveform(max_1_idx:min_idx)-halfmax));  half_1_idx = half_1_idx+max_1_idx-1; %adjust for range
    [~,half_2_idx] = min(abs(waveform(min_idx:max_2_idx)-halfmax));  half_2_idx = half_2_idx+min_idx-1;
    
%% Output variables
    ptIdx = [max_1_idx, min_idx, max_2_idx, half_1_idx, half_2_idx];
    ptVal = [max_1_val, min_val, max_2_val, halfmax]; 

end
