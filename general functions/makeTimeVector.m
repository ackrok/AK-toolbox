function timeVector = makeTimeVector(sigLength, Fs)
%% makeTimeVector
% Generate a vector with time in seconds matching the duration of
% photometry signal recording and based on sampling frequency of recording.
%
% timeVector = makeTimeVector(signalLength, Fs)
%
% INPUTS
% 'sigLength' - length of photometry signal (aka number of frames)
%       sigLength = length(FP);
%
% 'Fs' - sampling frequency, in Hz
%
% OUTPUTS
% 'timeVector' - vector with time, matching duration and sampling frequency
% of photometry recording

timeVector = [1:sigLength]./Fs; 
timeVector = timeVector(1:sigLength);

end