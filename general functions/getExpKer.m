function [expKer,kerSize] = getExpKer(tau_decay,tau_rise,Fs)
%getExpKer - Get exponential kernel
%
%   Usage:
%       [expKer,kerSize] = getExpKer(tau_rise,tau_decay,Fs);
%
%   Description: This function creates an exponential kernal that mimics
%   the dynamics (rise time and decay time) of a fluorescent sensor. The
%   function takes the dynamics and the sampling rate of spike times as the
%   input.
%
%   Input:
%       tau_rise - Rise time of sensor dynamics
%       tau_decay - Decay time of sensor dynamics
%       Fs - Sampling rate of spike times
%
%   Output:
%       expKer - Exponential Kernel
%       kerSize - The size of kernel in samples
%
%   Author: Pratik Mistry, 2020

tau_decay = (tau_decay/1000)*Fs; %Convert from ms into samples
expSize = tau_decay*5; % Approximate the size of the decay
kerSize = 2*expSize; % Specify the total size of the kernel
expKer = zeros(uint32(kerSize),1); 
expKer(uint32(expSize):end) = exp(-[0:expSize]/tau_decay); %Add the decay to the kernel
if tau_rise ~= 0 || ~isempty(tau_rise)
    tau_rise = (tau_rise/1000)*Fs; %Convert from ms into samples
    expKer(expSize-tau_rise:expSize) = ...
        linspace(0,1,length([expSize-tau_rise:expSize])); %Add the rise to the kernel
end

end