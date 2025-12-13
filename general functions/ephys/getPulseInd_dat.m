function pulseInd = getPulseInd_dat()
[pulseFile,fPath] = uigetfile('*.dat','Select Stimulation dat File');
fid = fopen(fullfile(fPath,pulseFile));
pulse = fread(fid,'int16');
pulseInd = getPulseOnsetOffset(pulse,max(pulse)/2);
end