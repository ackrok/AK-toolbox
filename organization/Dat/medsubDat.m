function medsubDat(nChan,memToUse)
%medsubData - Median Subtraction from dat file
%
%   Usage:
%       medsubDat(nChan,memToUse)
%
%   Description: This function takes dat files containing the raw data from
%   in vivo electrophysiology recordings and performs a median subtraction
%   on the data to decrease the effect of light artifacts in
%   opto-stimulation / opto-tagging experiments. The function takes
%   the number of channels in the file and the amount of memory to use as
%   inputs and then the function prompts the user to select the original
%   dat file and then input the new file name.
%
%   Input:
%       nChan - Number of channels in the dat file
%       memToUse - Amount of memory to use for the chunks being inputted
%   
%
%   Author: Pratik Mistry, 2020
%

[datName,datPath] = uigetfile('*.dat','Select Original Dat file');
datName_new = inputdlg('Enter new filename for dat file');
datName_new = datName_new{1};
[~,datNewExt] = strtok(datName_new,'.');
if (isempty(datNewExt))
    datName_new = [datName_new,'.dat'];
end
datFile = fullfile(datPath,datName);
datMap = getDatMemMap(datFile,nChan);
nSamp = datMap.Format{2}(2);
buff = (memToUse*(10^9))/(16*nChan);
nChunks = ceil(nSamp/buff);

fidout = fopen(fullfile(datPath,datName_new),'w');
fopen(fidout);

for ii = 1:nChunks
    fprintf('Subtracting median from chunk %d out of %d...\n',ii,nChunks);
    ind1 = int32(((ii-1)*buff)+1);
    if ii == nChunks
        ind2 = nSamp;
    else
        ind2 = int32(ii*buff);
    end
    tmpData = datMap.Data.dat(:,[ind1:ind2]);
    tmpData = reshape(tmpData,nChan,[]);
    tmpData = tmpData - median(tmpData);
    fwrite(fidout,tmpData,'int16');
end
fclose(fidout);
fprintf('Finished subtracting the median from data set!\n');
end