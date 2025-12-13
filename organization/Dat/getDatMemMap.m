function [datMemMap,datSize] = getDatMemMap(datFF,nChan)
%datMemMap - Create a memory map to dat files
%
%   Usage:
%       datMemMap = getDatMemMap(datFF,nChan);
%
%   Description: This function uses the full path to the dat file including
%   the number of channels and creates a memory map to it. This allows
%   users to access the dat file without having to read chunks using fread.
%   
%   Input:
%       datFF - Full path to the dat file
%       nChan - Number of channels in the dat file
%
%   Output:
%       datMemMap - Memory map data structure
%       datSize - A 1x2 array containing the size of data in the dat file
%
%   Author: Pratik Mistry, 2020
%

datStruct = dir(datFF);
nBytes = numel(typecast(cast(0,'int16'),'uint8'));
nSamp = datStruct.bytes/(nChan*nBytes);
datSize = [nChan,nSamp];
datMemMap = memmapfile(datFF,'Format',{'int16',[nChan nSamp],'dat'});
end