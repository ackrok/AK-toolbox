%% Generating cannula output structure
%
% Anya Krok, September 2022
%

%% Select .mat files you want to add to summary data structu
[fName,fPath] = uigetfile('*.mat','Select data files (multiple infusions) to add to cannula structure','MultiSelect','On');
if ~iscell(fName); fName = {fName}; end
fName = sort(fName);
beh = extractBeh(fPath, fName); % extract data

%% Organize data into cannula structure
infWindow = [10 40]; % window to analyze after infusion, in minutes
infOptions = {'saline','iGluR-antag','nAChR-antag','mAChR-antag','D1/2R-antag','D2R-ago','muscimol','other'};
for x = 1:length(beh)
    choice = menu(sprintf('Select infusion: %s',fName{x}),infOptions); % select data type
    beh(x).inf = infOptions{choice}; % store infusion label into structure
    if choice == 8
        infThis = inputdlg({'Enter infusion:'},'Input',1,{''});
        beh(x).inf = infThis{1};
    end
    
end
uni = unique({beh.inf}); % unique infusions
uni = fliplr(uni); % flip so that saline is first

cannula = struct; % initialize new structure
for y = 1:length(uni)
    cannula(y).s = beh(strcmp({beh.inf},uni{y}));
    cannula(y).inf = uni{y};
    cannula(y).win = infWindow.*ones(length(cannula(y).s),2); % window to analyze after infusion, in minutes
end

%% Save if indicated
choice = menu('Save cannula structure for future use?','Yes','No');
switch choice
    case 1
        if length(uni) > 1
            canName = inputdlg('Enter name for file to save cannula structure:','Input',1,{['achda_beh_cannulaDLS_',uni{1},'+',uni{2},'.mat']});
        else
            canName = inputdlg('Enter name for file to save cannula structure:','Input',1,{['achda_beh_cannulaDLS_',uni{1},'.mat']});
        end
        save(fullfile(fPath,canName{1}),'cannula');
end