function [] = createMoveEvt(mouserec,onSet,offSet,type)
%Generate .mov.evt file for Neuroscope
%Created By: Pratik Mistry
%Edited By: Anya Krok
%
% [] = createMoveEvt(mouserec,onSet,offSet,type)
%
% INPUT
%   'mouserec' - name of in vivo recording, to be used to name .mov.evt file
%   'onSet' - vector of event onset indices in *sec*
%   'offSet' - vector of event OFFset indices in *sec*
%   'type' - 'wheel' or 'opto'
%
% OUTPUT
%   .mov.evt file saved into 'R:\In Vivo\event-files\'
%

    path = ['R:\In Vivo\event-files\'];            %path = uigetdir([],'Select save location');

    switch type
        case 'wheel'
            file = cellstr([mouserec,'_wheel.mov.evt']);       %file = inputdlg('Input filename format "name.mov.evt"');
            fid = fopen(fullfile(path,file{1}),'w');
            for i = 1:length(onSet) %Onset and Offset should be same length
                fwrite(fid,sprintf('%f Movement_Onset\n',onSet(i)*1000)); %multiply by 1000 to convert onset/offset times from seconds to ms
                fwrite(fid,sprintf('%f Movement_Offset\n',offSet(i)*1000));
            end

        case 'opto'
            file = cellstr([mouserec,'_opto.mov.evt']);       %file = inputdlg('Input filename format "name.mov.evt"');
            fid = fopen(fullfile(path,file{1}),'w');
            for i = 1:length(onSet) %Onset and Offset should be same length
                fwrite(fid,sprintf('%f Pulse_Onset\n',onSet(i)*1000)); %multiply by 1000 to convert onset/offset times from seconds to ms
                fwrite(fid,sprintf('%f Pulse_Offset\n',offSet(i)*1000));
            end
    end
    
    fclose(fid);
    
end