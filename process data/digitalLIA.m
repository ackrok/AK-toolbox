function [sig] = digitalLIA(modSig,refSig,freq,Fs,lpCut,filtOrder)
%Digital Lock-In Amplifier Demodulation
%
%   [sig,filtStruct] = digitalLIA(modSig,refSig,Fs,lpCut,filtOrder)
%
%   Description: This function demodulates a signal using a phase sensitive
%   detection method. The function takes two inputs: modulated signal and
%   the recorded reference signal. If the input reference signal is in
%   phase with the modulated signal (as measured by a cross-correlation),
%   the function will perform a demodulation that is phase dependent. If
%   the user does not input an in phase reference signal, the function will
%   perform a phase independent demodulation that utilizes a reference
%   signal 90-degrees out of phase with the input reference signal; then
%   the function will demodulate using a quadrature demodulation method.
%
%   Multiplying two sinusoids:
%
%   A*sin(w_1*t+phi_1) .* B*sin(w_2*t+phi_2) =
%       0.5*A*B*cos((w_1-w_2)*t+(phi_1-phi_2)) -
%       cos((w_1+w_2)*t+(phi_1+phi_2));
%
%   Input:
%   - modSig - Measured modulated signal
%   - refSig - Measured driving signal used to drive the LED
%   - freq - Modulation Frequency
%   - Fs - Sampling Frequency
%   - lpCut - Low-Pass CutOff Frequency
%   - filtOrder - Filter Order
%
%   Output:
%   - sig - "Demodulated" signal
%
%
%   Author: Pratik Mistry, 2019

%Check to see if signals are column vectors; if not, it corrects the
%signal orientation to column vectors
modSig = sub_chkSigSize(modSig); refSig = sub_chkSigSize(refSig);
%Check to see if the signals are the same length because the signal
%multiplication requires the signals to be exactly the same length
if (size(modSig,1) ~= size(refSig,1))
    disp('ERROR: Signals are not the same length. Signals need to be the same length for the code to perform phase sensitive detection');
    sig = 0;
    return;
else
    %Check to see if user inputted an in-phase reference signal
    %Perform a cross-correlation. If maximum correlation exists at 0,
    %then the signals are in phase.
    [x,lag] = xcov(modSig,refSig,'coeff');
    phaseDiff = lag(find(x==max(x)));
    %Clear variables to free space
    clear x lag;
    
    %Normalize Reference Signal and ensure the amplitude goes from +2
    %to -2V --> This step ensures that you are maintaining the original
    %ampltiude from the modulated photometry signal
    refSig = refSig-min(refSig); refSig = refSig/max(refSig); 
    refSig = refSig*4;
    
    bpFilt = designfilt('bandpassiir', 'FilterOrder', 6,...
        'HalfPowerFrequency1',freq-10, 'HalfPowerFrequency2',freq+10,...
        'SampleRate', Fs, 'Designmethod', 'butter');
    modSig = filtfilt(bpFilt, modSig);
   
    lpFilt = designfilt('lowpassiir','FilterOrder',filtOrder,...
        'HalfPowerFrequency',lpCut,'SampleRate',Fs,...
        'DesignMethod','butter');
    
    modSig = modSig - mean(modSig);
    refSig = refSig - mean(refSig);
    
    if phaseDiff == 0 %Signals are in-phase perform standard PSD
        PSD = modSig.*refSig;
        sig = filtfilt(lpFilt,PSD); % sig = PSD;
    else %Signals are not in-phase compute a quadrature using reference signal shift 90 degrees
        refSig_90 = [diff(refSig);refSig(end)-refSig(end-1)];
        PSD_1 = modSig.*refSig;
        PSD_1 = filtfilt(lpFilt,PSD_1);
        PSD_2 = modSig.*refSig_90;
        PSD_2 = filtfilt(lpFilt,PSD_2);
        sig = hypot(PSD_1,PSD_2);
    end
end

end

function adjSig = sub_chkSigSize(orgSig)
%Check Signal Size/Orientation
%
%   [adjSig] = sub_chkSigSize(orgSig);
%
%   Description: This sub-function ensures the orientation of the signals
%   are both column vectors. This function is necessary because if the
%   vectors are not the same orientation for the phase sensitive detection,
%   it will throw an error.
%
%   Input:
%   - orgSig - Original Signal
%
%   Output:
%   - adjSig - Adjusted Signal
%
%

%This function uses the size function to ensure that size(sig,2) = 1
%The second index in size is the number of columns
if size(orgSig,2) ~= 1
    adjSig = orgSig';
else
    adjSig = orgSig;
end
end