function time = makeTime(nSamp, Fs)
% Description: generate time vector using length of signal and sampling
% frequency
%
% time = makeTime(nSamp, Fs)
%
% INPUTS
% 'nSamp': value, number of samples in signal (eg, length(fp))
% 'Fs'   : value, sampling frequency in Hz
%
% OUTPUTS
% 'time' : column vector with time from 0 to length of recording, in sec
%
% Anya Krok, January 2026

time = linspace(0, nSamp/Fs, nSamp); % time vector
time = time(:); % make into column vector


end