function [norm, f, flog, p1_mat, rawS] = getCannula_fft(window, varargin)
%FFT to analyze the frequency components present in photometry signal for 
% cannula infusion experiments
%
%   [norm, f, flog, p1_mat, rawS] = getCannula_fft()
%   [norm, f, flog, p1_mat, rawS] = getCannula_fft(rawS)
%
% INPUTS
%   'rawS' OR if raw signals are not already loaded into workspace in a 
%   structure called rawS, user will then be prompted to select the .mat 
%   files that contain photometry data to be analyzed.
%
% OUTPUTS
%   'norm' - matrix with FFT output for all recordings after normalization
%   across range of 0.1 to 100 Hz
%   'f' - frequency domain vector
%   'flog' - log-scale frequency domain vector
%   'p1_mat' - matrix with FFT output
%   'rawS' - structure with raw photoemtry signals
%
% Anya Krok, July 2022

%% LOAD RAW SIGNALS INTO WORKSPACE
switch nargin
    case 2
        rawS = varargin{1};
        rawS = extractRaw_fft_cannula(window, rawS);
    case 1
        %% Select .mat files you want to add to summary data structu
        [fName,fPath] = uigetfile('*.mat','MultiSelect','On');
        rawS = extractRaw_fft_cannula(window, fPath, fName);
end
if ~exist('rawS','var'); error('No variable called rawS exists'); end

%% FFT
[p1_mat, f] = getFft(rawS);
fprintf('FFT done! \n');

%% Normalize FFT
range_norm = [0.1 100]; % Range for normalization of FFT output, in Hz
tmpNorm = [];
r = [find(f == range_norm(1)):find(f == range_norm(2))]; % Restrict to [0.01 100]
flog = log10(f(r)); % Log-scale frequency vector
for x = 1:size(p1_mat,2)
    a = log10(p1_mat(r,x));
    a_end = nanmean(log10(p1_mat(find(f == 40):find(f == 50),x)));
    vec_norm = (a - a_end)./(a(1) - a_end); 
    % vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100], scaling so range covers [0 1]
    tmpNorm(:,x) = vec_norm;
end
norm = tmpNorm;
fprintf('Normalization done! \n');

end