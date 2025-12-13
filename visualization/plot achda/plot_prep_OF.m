% pull = xlsread('C:\Users\Anya\Desktop\FP_LOCAL\openField\SummaryNov22_1_30Hz.xlsx');
% pull_2 = xlsread('C:\Users\Anya\Desktop\FP_LOCAL\openField\SummaryNov22_2_60Hz.xlsx')
load('C:\Users\Anya\Desktop\FP_LOCAL\openField\beh_OF_SummaryNov22.mat')

%% Check mov/imm
x = 1;
mov = zeros(length(beh(x).time),1); res = zeros(length(beh(x).time),1);
for y = 1:length(beh(x).on); mov(beh(x).on(y):beh(x).off(y)) = 1; end
for y = 1:length(beh(x).onRest); res(beh(x).onRest(y):beh(x).offRest(y)) = 1; end
figure; hold on
area(beh(x).time, mov, 'FaceColor', 'g', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.1);
area(beh(x).time, res, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.1);

%% Refine mov/imm to be >4s long
thres = 4; % Remove bouts that are less than X seconds long
for x = 1:length(beh)
    Fs = beh(x).Fs;
    len = (beh(x).off - beh(x).on)./Fs; % Length of bouts, in seconds
    rmv = find(len < thres); % Index of bouts that are less than X seconds long
    beh(x).on(rmv) = []; beh(x).off(rmv) = []; % Remove short bouts
    len = (beh(x).offRest - beh(x).onRest)./Fs; % Length of bouts, in seconds
    rmv = find(len < thres); % Index of bouts that are less than X seconds long
    beh(x).onRest(rmv) = []; beh(x).offRest(rmv) = []; % Remove short bouts
end

%% AK231 232 233 234 236
sigEdge = 15; % Photometry analysis removes first and last 15s of recording
%{
for x = 1:5
    cam = beh(x).cam(:);
    Fs = length(cam)/1800; % Sampling freq of camera
    
    a = 2*x;
    tmp = pull(:,a); % Extract from matrix
    tmp(length(tmp)-Fs*sigEdge+1:end) = []; % Remove last 15s
    tmp(1:Fs*sigEdge) = []; % Remove first 15s, now length matches photometry
    [start, len, k1] = ZeroOnesCount(tmp); % Start and length of sequence of 1's
    start = start(:); len = len(:);
    mov = [];
    mov(:,1) = start(1:k1);
    mov(:,2) = mov(:,1) + len(1:k1) - 1;
    
    a = 2*x + 1;
    tmp = pull(:,a); % Extract from matrix
    tmp(length(tmp)-Fs*sigEdge+1:end) = []; % Remove last 15s
    tmp(1:Fs*sigEdge) = []; % Remove first 15s, now length matches photometry
    [start, len, k1] = ZeroOnesCount(tmp); % Start and length of sequence of 1's
    start = start(:); len = len(:);
    res = [];
    res(:,1) = start(1:k1);
    res(:,2) = res(:,1) + len(1:k1) - 1;
    
    mov = cam(mov); % Adjust to be in samples rather than camera time stamps
    res = cam(res); % Adjust to be in samples rather than camera time stamps 
    if mov(1,1) == 0; mov(1,:) = []; end
    if res(1,1) == 0; res(1,:) = []; end

    beh(x).on = (mov(:,1)); 
    beh(x).off = (mov(:,2));
    beh(x).onRest = (res(:,1)); 
    beh(x).offRest = (res(:,2));
end
%}

%% AK237 238
sigEdge = 15; % Photometry analysis removes first and last 15s of recording
%{
for x = 1:2
    cam = beh(x+5).cam(:);
    Fs = length(cam)/1800; % Sampling freq of camera
    
    a = 2*x;
    tmp = pull_2(:,a); % Extract from matrix
    tmp(length(tmp)-Fs*sigEdge+1:end) = []; % Remove last 15s
    tmp(1:Fs*sigEdge) = []; % Remove first 15s, now length matches photometry
    [start, len, k1] = ZeroOnesCount(tmp); % Start and length of sequence of 1's
    start = start(:); len = len(:);
    mov = [];
    mov(:,1) = start(1:k1);
    mov(:,2) = mov(:,1) + len(1:k1) - 1;
    
    a = 2*x + 1;
    tmp = pull_2(:,a); % Extract from matrix
    tmp(length(tmp)-Fs*sigEdge+1:end) = []; % Remove last 15s
    tmp(1:Fs*sigEdge) = []; % Remove first 15s, now length matches photometry
    [start, len, k1] = ZeroOnesCount(tmp); % Start and length of sequence of 1's
    start = start(:); len = len(:);
    res = [];
    res(:,1) = start(1:k1);
    res(:,2) = res(:,1) + len(1:k1) - 1;
    
    mov = cam(mov); % Adjust to be in samples rather than camera time stamps
    res = cam(res); % Adjust to be in samples rather than camera time stamps 
    if mov(1,1) == 0; mov(1,:) = []; end
    if res(1,1) == 0; res(1,:) = []; end

    beh(x+5).on = (mov(:,1)); 
    beh(x+5).off = (mov(:,2));
    beh(x+5).onRest = (res(:,1)); 
    beh(x+5).offRest = (res(:,2));
end
%}