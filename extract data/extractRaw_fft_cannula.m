function [rawS] = extractRaw_fft_cannula(window, varargin)
%%Extract raw photometry data from multiple recordings into a single
%%structure and parse by specified behavioral state, for FFT analysis with
%%function getFft
%
% [rawS] = extractRaw_fft_cannula(window) - extract new data into structure
% [rawS] = extractRaw_fft_cannula(window, fPath, fName) - extract new data into structure
% [rawS] = extractRaw_fft_cannula(window, rawS) - only parse already extracted data
%
% Description: Extract raw photometry signal from multiple recording files 
% into a larger structure to be used for FFT analysis with function getFft
% User will be prompted to select photometry signal and behavioral state 
% to be used. Recordings lacking behavioral state (e.g. no immobility) 
% will not be included in the output structure.
% Cannula: restrict signal to post-infusion window, default +20 to +40min
%
% INPUTS
%   'window' - Window, in seconds, used for analyzing infusion effect
%   'fPath' - Character array containing folder path where data files are
%       example: 'R:\tritsn01labspace\Anya\FiberPhotometry\AK201-206\220105'
%   'fName' - Cell array, with each cell containing file names for each
%               recording to be added to structure
%   'rawS' - Structure previously generated using extractRaw_fft function
%
% OUPUTS
%   'rawS' - Structure with raw photometry signals from multiple recordings
%       rawS(x).fp_sub, field containing parsed photometry signal during
%       specified behavioral state
%
% Anya Krok, January 2022
%
    %% INPUTS
    switch nargin
        case 1
            [fName,fPath] = uigetfile('*.mat','MultiSelect','On');
            if ~iscell(fName); fName = {fName}; end
            fName = sort(fName);
        case 3
            fPath = varargin{1};
            fName = varargin{2};
            if ~iscell(fName); fName = {fName}; end
            fName = sort(fName);
        case 2
            rawS = varargin{1};
    end
    behState = menu('Select behavioral state','Immobility','Locomotion','Reward','Full Trace');
    y = menu('Select photometry signal','ACh (1)','DA (2)');
    %% Extract data
    if nargin ~= 2 % If need to load data into new structure
        rawS = struct;
        h = waitbar(0, 'Extracting raw photometry signals into structure');
        for f = 1:length(fName)
            fprintf('Extracting raw photometry data %s ... ',fName{f});
            load(fullfile(fPath,fName{f})); % Load raw data file
            [an,b] = strtok(fName{f},'_'); [day,b] = strtok(b,'_'); [inf,~] = strtok(b,'_'); % Parse file name
            x = 1 + length(rawS);
            rawS(x).rec = [an,'-',day]; 
            rawS(x).site = 'DLS'; % CHANGE or remove
            rawS(x).inf = inf;

            %% Pull parameters required for this analysis
            if isfield(data.gen,'params')
                params = data.gen.params; % Extract params structure
                dsRate = params.dsRate; 
                dsType = params.dsType; % General downsampling parameter
                rawFs = data.gen.acqFs; 
                Fs = data.gen.Fs;
            else
                error('No parameters saved during processData');
            end

            %% Extract photometry and behavior data
            rawS(x).FPnames = data.acq.FPnames;
            rawS(x).rawFP = data.acq.FP;
            rawS(x).rawFs = rawFs;
            if isfield(data,'final')
                if isfield(data.final,'mov')
                    if ~isempty(data.final.mov.onsets)
                        rawS(x).on = data.final.mov.onsets.*dsRate;
                        rawS(x).off = data.final.mov.offsets.*dsRate;
                    end
                    if ~isempty(data.final.mov.onsetsRest)
                        rawS(x).onRest = data.final.mov.onsetsRest.*dsRate;
                        rawS(x).offRest = data.final.mov.offsetsRest.*dsRate;
                    end
                end
            end
            if isfield(data.final,'rew')
                if isfield(data.final.rew,'onset')
                    rawS(x).reward = data.final.rew.onset.*dsRate;
                end
            end
            win_inf = window.*rawFs; % convert infusion window into samples
            fprintf('DONE.\n');
            waitbar(f/length(fName),h);

        end
        close(h);
        if isempty(rawS(1).rawFs); rawS(1) = []; end
    end

    %% Parse photometry signal during specified behavioral state
%     behState = menu('Select behavioral state','Immobility','Locomotion','Reward','Full Trace');
%     y = menu('Select photometry signal',rawS(1).FPnames);
    rmv = zeros(length(rawS),1);
    switch behState
        case 1
            fprintf('Extracting signal (this will take a while!) ... \n')
            h = waitbar(0,'Extracting signal during behavioral state');
            for x = 1:length(rawS)
                rawFP = rawS(x).rawFP{y}; rawFs = rawS(x).rawFs; % photometry signal
                nSampRaw = length(rawS(x).rawFP{y}); % number of samples
                win_inf = window.*rawFs; % convert infusion window into samples
                if isfield(rawS,'reward')
                    rewWindow = rawFs;
                    rew = rawS(x).reward(rawS(x).reward > win_inf(1) & rawS(x).reward < win_inf(2));
                    idx_rew = extractEventST([1:nSampRaw]', floor(rew), floor(rew)+rewWindow, 1); % identify recording indices during reward
                else; idx_rew = [];
                end
                onRest = rawS(x).onRest(rawS(x).onRest > win_inf(1) & rawS(x).onRest < win_inf(2));
                offRest = rawS(x).offRest(rawS(x).onRest > win_inf(1) & rawS(x).onRest < win_inf(2));
                idx_imm = extractEventST([1:nSampRaw]', onRest, offRest, 1); % identify recording indices during immobility
                idx_imm = idx_imm(idx_imm > win_inf(1) & idx_imm < win_inf(2));
                if isempty(idx_imm)
                    fprintf('%s - no immobility \n', rawS(x).rec); 
                    rmv(x) = 1; % change index to be 1, this recording will be removed
                end
                idx_imm = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward from immobility indices
                rawS(x).fp_sub = rawFP(idx_imm); % extract signal during immobility
                rawS(x).fp_lbl = rawS(x).FPnames{y};
                rawS(x).behState = 'immobility';
                waitbar(x/length(rawS),h);
            end; close(h);
%         case 2
%             fprintf('Extracting signal during %s (this will take a while!) ...',behState)
%             h = waitbar(0,'Extracting signal during behavioral state');
%             for x = 1:length(rawS)
%                 nSampRaw = length(rawS(x).rawFP{y});
%                 if isfield(rawS,'reward')
%                     rewWindow = rawS(x).rawFs;
%                     idx_rew = extractEventST([1:nSampRaw]', floor(rawS(x).reward), floor(rawS(x).reward)+rewWindow, 1); % identify recording indices during reward
%                 else; idx_rew = [];
%                 end
%                 idx_loc = extractEventST([1:nSampRaw]', rawS(x).on, rawS(x).off, 1); % identify recording indices during locomotion
%                 if isempty(idx_loc)
%                     fprintf('%s - no locomotion \n', rawS(x).rec); 
%                     rmv(x) = 1; % change index to be 1, this recording will be removed 
%                 end
%                 idx_loc = idx_loc(~ismember(idx_loc, idx_rew)); % exclude reward from locomotion indices
%                 rawS(x).fp_sub = rawS(x).rawFP{y}(idx_loc); % extract signal during locomotion
%                 rawS(x).behState = 'locomotion';
%                 waitbar(x/length(rawS),h);
%             end; close(h);
%         case 3
%             fprintf('Extracting signal during %s (this will take a while!) ...',behState)
%             h = waitbar(0,'Extracting signal during behavioral state');
%             for x = 1:length(rawS)
%                 nSampRaw = length(rawS(x).rawFP{y});
%                 if isfield(rawS,'reward')
%                     rewWindow = rawS(x).rawFs;
%                     idx_rew = extractEventST([1:nSampRaw]', floor(rawS(x).reward), floor(rawS(x).reward)+rewWindow, 1); % identify recording indices during reward
%                     rawS(x).fp_sub = rawS(x).rawFP{y}(idx_rew); % extract signal during reward
%                     rawS(x).behState = 'reward';
%                 else
%                     error('no reward');
%                 end
%                 waitbar(x/length(rawS),h);
%             end; close(h);
        case 4
            for x = 1:length(rawS)
                rawS(x).fp_sub = rawS(x).rawFP{y}; % full trace
                rawS(x).fp_lbl = rawS(x).FPnames{y};
                rawS(x).behState = 'full';
            end
    end
    rawS(rmv == 1) = []; % remove recordings where no fp_sub extracted
    
end