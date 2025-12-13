%% Convert Scan Image Files
clear all

data = struct('mouse',[],'date',[],'acq',struct());
data.acq = struct('FPnames','','nFPchan',[],'FP',[]);

%% Pull Wheel Data
choice = menu('Are you adding Wheel Data?','Yes','No');

switch choice
    case 1
        [fFiles,fPath] = uigetfile('*.mat','MultiSelect','On');
        if (~iscell(fFiles))
            fFiles = {fFiles};
        end
        fFiles = sort(fFiles);
        for n = 1:length(fFiles)
            [fName,~] = strtok(fFiles{n},'.');
            wheel = getfield(open([fPath,fFiles{n}]),fName);
            data.acq(n).wheel = wheel.data;
        end
    case 0 || 2
        disp('No Wheel Data Entered');
end

%% Pull Photometry Data
choice = menu('Are you adding Photometry Data?','Yes','No');
switch choice
    case 1
        nFP = 1;
        FPname = inputdlg('Enter Name of FP Channel');
        [fFiles,fPath] = uigetfile('*.mat','MultiSelect','On');
        if (~iscell(fFiles))
            fFiles = {fFiles};
        end
        fFiles = sort(fFiles);
        for n = 1:length(fFiles)
            [fName,~] = strtok(fFiles{n},'.');
            FP = getfield(open([fPath,fFiles{n}]),fName);
            data.acq(n).FP{nFP,1} = FP.data;
            data.acq(n).FPnames{nFP,1} = FPname;
        end
        while(1)
            choice_2 = menu('Do you want to add another Photometry channel?','Yes','No');
            switch choice_2
                case 1
                    nFP = nFP+1;
                    [fFiles,fPath] = uigetfile('*.mat','MultiSelect','On');
                    if (~iscell(fFiles))
                        fFiles = {fFiles};
                    end
                    for n = 1:length(fFiles)
                        [fName,~] = strtok(fFiles{n},'.');
                        FP = getfield(open([fPath,fFiles{n}]),fName);
                        data.acq(n).FP(nFP,1) = FP.data;
                    end
                case 0 || 2
                    break;
            end
        end
    case 0 || 2
        disp('No Photometry Data to Enter');
end

%% Pull Opto Pulse Data
choice = menu('Are you adding Opto Pulses?','Yes','No');

switch choice
    case 1
        [fFiles,fPath] = uigetfile('*.mat','MultiSelect','On');
        if (~iscell(fFiles))
            fFiles = {fFiles};
        end
        fFiles = sort(fFiles);
        for n = 1:length(fFiles)
            [fName,~] = strtok(fFiles{n},'.');
            opto = getfield(open([fPath,fFiles{n}]),fName);
            data.acq(n).opto = {opto.data};
        end
    case 0 || 2
        disp('No Opto Pulses Entered');
end

%% Enter General Info
[genInfo] = inputdlg({'Enter Mouse Name','Enter Acquisition Date','Enter Sampling Frequency'});
data.mouse = genInfo{1};
data.date = genInfo{2};
for n = 1:length(data.acq);
    data.acq(n).nFPChan = length(data.acq(n).FP);
    data.acq(n).Fs = str2double(genInfo{3});
end

clearvars -except data fPath

newName = inputdlg('Enter New File Name');

save([fPath,newName{1},'.mat'],'data');