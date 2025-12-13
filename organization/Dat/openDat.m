function [datSignal, fName, fPath] = openDat(varargin)
%Extract signal from single channel .dat file
%
% [datSignal, fName, fPath] = openDat()
% [datSignal, fName, fPath] = openDat(datID)
%
% INPUT
%   'datID' (optional) - string with file ID, which will be displayed 
%       in file selection prompt 
%
% OUTPUT
%   'datSignal' - vector with signal from .dat file
%   'fName' - .dat file name
%   'fPath' - .dat file path
%
% Created by: Pratik Mistry, edited by Anya Krok

switch nargin
    case 0
        [fName, fPath] = uigetfile('*.dat','Select .dat file to extract');
    
    case 1
        datID = varargin{1};
        [fName, fPath] = uigetfile('*.dat',['Select .dat file to extract: ',datID]);
        
    case 2
        fPath = varargin{1}; fName = varargin{2};
end

datFF = [fPath, fName]; %dat full file name
fid = fopen(datFF); 
datSignal = fread(fid,'int16'); %read from file identifier (fid), must specify 'int16'
fclose(fid);

end