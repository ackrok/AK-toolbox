function cost = van_rossum(st1,st2,tau,Fs)
%van_rossum - Get van rossum spike distance cost
%
%   Usage:
%       [cost] = van_rossum(st1,st2,tau,Fs)
%
%   Description: This function takes two spike trains as an input, 
%   convolves them with an exponential kernel to obtain a continuous
%   representation of the spike trains. These new continuous signals are
%   used to calculate a 'cost' metric to quantify similarity between two 
%   spike trains
%
%   Cost = sqrt(1/tau * sum((x - y)^2)
%       where x and y are the convolved spike times
%
%   Input:
%       st1 - Spike train 1
%       st2 - Spike train 2
%       tau - Kernel tau
%       Fs - Sampling rate
%
%   Output:
%       cost - Metric for spike train similarity
%
%   Author: Pratik Mistry, 2020

    [exp_ker, ~] = getExpKer(tau,0,Fs);
    st1 = uint32(st1*Fs); st2 = uint32(st2*Fs);
    if max(st1) > max(st2)
        Ls = max(st1);
    else
        Ls = max(st2);
    end
    
    st1_conv = zeros(Ls,1); st2_conv = zeros(Ls,1);
    st1_conv(st1) = 1; st2_conv(st2) = 1;
    %st1_conv = conv(st1_conv,exp_ker); st2_conv = conv(st2_conv,exp_ker);
    st1_conv = nanconvn(st1_conv,exp_ker); st2_conv = nanconvn(st2_conv,exp_ker);
    
    cost = (st1_conv - st2_conv).^2;
    cost = sqrt(sum(cost) / (tau/1000));
end