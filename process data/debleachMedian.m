function [sig_norm, sig_iter] = debleachMedian(signal, Fs)
% Debleaching photometry signal with two methods
%
% Description: perform local detrending by subtracting local average from
% staggered blocks throughout signal. First pass uses non-overlapping
% blocks and mean-centers trace - output is 'sig_norm'. Second pass runs
% multiple staggered passes that start at shifted offsets and
% median-centers to combine, thus produces overlapping coverage across
% session to reduce block-edge bias - output is 'sig_iter'.
% Note: analysis is session-centric not trial-based. 
%
% [sig_norm, sig_iter] = debleachMedian(signal, Fs)
%
% INPUTS
%   'signal': column vector with raw photometry signal
%   'Fs': sampling frequency in Hz
%
% OUTPUTS
%   'sig_norm': mean-centered non-overlapping blocks normalized by stdev
%   'sig_iter': median-centered staggered blocks normalized by median
%               absolute deviation
%
% Anya Krok, March 2026
% Adapted from Bruno et al. 2021 - pMAT fiber photometry analysis FP_DEBLEACHED

%% 
windowSec = 0.5; % window duration in seconds for local centering
iterations = 3; % number of staggered iterations for computing baseline

N = numel(signal); 
blockLen = max(1, round(windowSec * Fs)); % block length, in samples
% incremental shift between staggered iterations:
cut = max(1, round(blockLen / iterations)); 

%%
% First pass steps through non-overlapping blocks. For each block, compute
% local mean and median. Use local mean to calculate mean-centered trace
sig_centered = nan(N, 1);
for idx = 1 : blockLen : N
    i2 = min(N, idx + blockLen - 1); % end index of current block
    block = signal(idx:i2); % extract block segment
    mu = mean(block, 'omitnan'); % mean
    sig_centered(idx:i2) = signal(idx:i2) - mu;
end
sigma = std(sig_centered, 0, 'omitnan'); % std across entire mean-centered trace
if sigma == 0, sigma = 1; end % avoid dividing by zero
sig_norm = sig_centered ./ sigma; % normalize mean-centered trace by stdev

%%
% Second pass is to compute median-centered blocks with staggered
% iterations. Each pass substracts local medians for each block that 
% starts at a different offset. This will reduce sensitivity to block edges.
iterAll = nan(N, iterations-1);

for it = 1:(iterations-1)
    iterThis = nan(N,1); % preallocate current iteration's median-centered result
    startIdx = (it * cut); % starting index for this iteration
    % iterate blocks starting at startIdx, then every blockLen samples.
    % this will create overlapping coverage across iterations
    for idx = startIdx:blockLen:N
        if idx <= 0
            continue;
        end
        i1 = idx;
        i2 = min(N, idx + blockLen - 1);
        if i1 > N % If iteration starts beyond N, break
            break;
        end
        block = signal(i1:i2);
        med = median(block, 'omitnan'); % local median
        iterThis(i1:i2) = signal(i1:i2) - med; % median adjustment
    end
    % handle leading region (if startIdx > 1) when staggered start does not
    % cover initial samples -- for samples before startIdx will fill them
    % with median of earliest segment to avoid leaving leading NaNs 
    if startIdx > 1
        leadEnd = min(N, startIdx-1);
        medLead = median(signal(1:leadEnd), 'omitnan'); % median of 1st segment
        iterThis(1:leadEnd) = signal(1:leadEnd) - medLead;
    end
    iterAll(:, it) = iterThis; % store in output matrix
end

% Lastly, combine staggered iterations and normalize by median.
% use sample-wise median across iterations to reduce block-edge artifacts
med_tmp = median(iterAll, 2, 'omitnan'); % median across columns for each sample
% then normalize by MAD (median absolute deviation) to create z-measure
mad_tmp = mad(med_tmp, 1); 
if mad_tmp == 0; mad_tmp = eps; end % avoid dividing by zero (eps is tiny positive)
sig_iter = med_tmp ./ mad_tmp; 

end
