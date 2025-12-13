function [ref] = findRef(modFreq,refSig,Fs)
%Find reference signal
%
%   [ref] = findRef(modFreq,refSig,Fs)
%
%   Description: This function goes through all the reference signals in
%   the refSig cell and tries to find the signal that corresponds to the
%   input modulation frequency using an FFT with a +/-5Hz window
%
%   Input:
%   - modFreq - Modulation frequency specified by the user
%   - refSig - Cell array containing reference signals
%   - Fs - Sampling frequency of the traces
%
%   Output:
%   - ref - Reference signal that has the same frequency as modFreq
%
%   Author: Pratik Mistry 2019

for n = 1:length(refSig)
    tmpRef = refSig{n}; %Pull the signal
    [refMag,refFreq] = calcFFT(tmpRef-mean(tmpRef),Fs); %Calculate the frequency and find the x and y axis
    maxRefFreq = refFreq(find(refMag == max(refMag))); %Find the frequncy that the max power occurs
    %This if-statement checks to see if maxRefFreq matches the inputted
    %modulation frequency
    if modFreq >= (maxRefFreq-5) && modFreq <= (maxRefFreq+5)
        ref = tmpRef;
    end
end

end

function [fftamp,fftfreq,fftphase] = calcFFT(y,fs,varargin)
%calcFFT - Calculate the Fourier Power Spectrums
%   Created By: Pratik Mistry
%   Created On: 31 January 2019
%   Edited On: 5 June 2019
%
%   [fftamp,fftfreq,fftphase] = calcFFT(y,fs)
%
%   Description: This function will calculate the Fourier spectrum of a
%   signal and return the phase and magnitude spectrums. The function can
%   also plot the variables if necessary
%
%   Input:
%   - y - Input signal
%   - fs - Sampling frequency
%
%   Output:
%   - fftamp - Vector of amplitudes for positive frequencies
%   - fftphase - Phase associated with positive frequencies
%   - fftfreq - Vector of positive freqencies -- Xaxis
%
% Author: Pratik Mistry, 2019

Y=fft(y); %Calculate the the Fast-Fourier Transform of the desired trace
L=length(Y); %Get the length of the calculated transform
Y_pos=Y(1:(L/2)); %Since Fourier spectrograms are symmetrical, we only take the first half of the transformed signal
fftfreq=(0:(L/2)-1)*fs/L; %Obtain the frequency vector used for plotting
fftamp=abs(Y_pos); %Since FFT's are complex in nature, we only want to take the magnitude of the values
Y_pos(abs(Y_pos)<1e-6) = 0; %Remove any super low power contributions

end