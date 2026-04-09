function [p, rawP, adjP, isParametric] = checkParametric(mat)
%% checkParametric
% Before choosing parametric vs non-parametric testing,
% (a) Inspect residuals
%   use Shapiro–Wilk or Lilliefors for normality of differences.
% (b) Test sphericity (Mauchly) if >2 repeated levels; 
%   if violated, apply Greenhouse–Geisser/Hyunh–Feldt corrections 
%   (fitrm/ranova supports epsilons).%
%
% Syntax:
%   [p, rawP, adjP, isParametric] = checkParametric(mat);
%
% Input:
%   mat - matrix, n x k (n animals, k conditions)
%
% Output:
%   p    - p-value from repeated-measures ANOVA or Friedman
%   rawP - after paired comparisons 
%   adjP - after Bonferroni correction
%   isParamteric - true if normality and sphericity holds
%
% Decision Flow:
% If pairwise differences ~ normal AND sphericity holds 
%   → use repeated‑measures ANOVA.
% If differences ~ normal but sphericity violated 
%   → use RM‑ANOVA but report GG corrected p.
% If differences nonnormal (especially across several pairs) 
%   → use Friedman's test.
%
% Anya Krok, March 2026

%%
n = size(mat,1);
k = size(mat,2);

%% Inspect residuals / normality of differences
% For k=3, check normality of the differences for each pair. 
% Use Lilliefors (lillietest) 
% or Shapiro–Wilk (if you have an implementation). 
% Also inspect QQ plots and histograms.
fprintf('\n \n Check parametric or non-parametric. \n \n');
fprintf('(a)');
% M: n x k matrix
pairs = nchoosek(1:k, 2);
alpha = 0.05;
normality = [];
for a = 1:size(pairs,1)
    i = pairs(a,1); j = pairs(a,2);
    d = mat(:,i) - mat(:,j);
    % Visual checks
    % figure;
    % subplot(1,2,1); histogram(d); title(sprintf('Diff %d-%d', i, j));
    % subplot(1,2,2); qqplot(d); title(sprintf('QQ %d-%d', i, j));
    % Normality test: Lilliefors (Kolmogorov-Smirnov with estimated params)
    [h,pval] = lillietest(d, 'Alpha', alpha);
    fprintf('Pair %d-%d: Lilliefors h=%d, p=%.3g\n', i, j, h, pval);
    % If you have a Shapiro-Wilk function (swtest), prefer it:
    % [h_sw,p_sw] = swtest(d, alpha);
    normality(a) = h;
end

% Interpret: h==0 means fail to reject normality (OK). 
% If several differences violate normality, 
% consider nonparametric approach (Friedman + pairwise signrank).

%% Repeated measures ANOVA and sphericity (Mauchly)
% Use fitrm + ranova and run mauchly to test sphericity. 
% ranova returns p-values with Greenhouse–Geisser (pValueGG) 
% and Huynh–Feldt (pValueHF) corrections.

% Build table with variable names for conditions
T = array2table(mat, 'VariableNames', {'Cond1','Cond2','Cond3'});
% Within design: indicate measurement levels (1..k)
Within = table([1;2;3], 'VariableNames', {'Time'});  % simple within factor
rm = fitrm(T, 'Cond1-Cond3 ~ 1', 'WithinDesign', Within);

% Mauchly's test
mTbl = mauchly(rm);
sphericity = mTbl.pValue; 
fprintf('(b) Mauchly test: p-value = %1.3f \n', sphericity);
% disp(mTbl);   % contains pValue for sphericity test

% Interpretation:
% If mauchly pValue > alpha → sphericity holds; use ranova pValue.
% If mauchly pValue ≤ alpha → sphericity violated; use pValueGG or pValueHF 
%   from ranova (Greenhouse–Geisser is conservative; 
%   Huynh–Feldt is less conservative).

%% decision flow
% If pairwise differences ~ normal AND sphericity holds 
% → use repeated‑measures ANOVA.
% If differences ~ normal but sphericity violated 
% → use RM‑ANOVA but report GG/HF corrected p (ranova provides them).
% If differences nonnormal (especially across several pairs) → use Friedman's test
%
% Repeated measures ANOVA with epsilon-corrected p-values
if ~all(normality)
    ranovatbl = ranova(rm);
    % disp(ranovatbl);  % contains pValue, pValueGG, pValueHF, pValueLB
    if sphericity > alpha
        p = ranovatbl.pValue(1);
    else
        p = ranovatbl.pValueGG(1);
        % Greenhouse–Geisser (GG) more conservative (reduces type 1 error)
        % better for smaller sample sizes compared to HF
    end
    fprintf('\n --> PARAMETRIC. \n \n');
    fprintf('Repeated-measures ANOVA: \n   p-value = %1.4f \n',p);

    if p < alpha 
    % If RM‑ANOVA significant, do paired comparisons with correction:
    %   Paired t-tests with Bonferroni/Holm, or
    %   use multcompare on the repeated model.
        % Pairwise paired t-tests and Bonferroni correction
        fprintf('\nPaired comparisons with correction: \n');
        rawP = zeros(k,1);
        for a = 1:size(pairs,1)
            [~, rawP(a)] = ttest(mat(:,pairs(a,1)), mat(:,pairs(a,2))); % parametric
            % rawP(a) = signrank(mat(:,pairs(a,1)), mat(:,pairs(a,2))); % nonparametric
        end
        adjP = min(rawP * numel(rawP), 1);
        fprintf('   Raw p: %s\n   Adj p (Bonferroni): %s \n \n', ...
            mat2str(round(rawP,3)'), mat2str(round(adjP,3)'));
        
        % multcompare requires a repeated factor name; for simple WithinDesign use:
        fprintf('Within-subject comparison: \n');
        results = multcompare(rm, 'Time');   % examines levels of the within-subject factor
        % disp(results);
    end

    % non-parametric
elseif all(normality) || sphericity <= alpha
    % If normality or sphericity is violated → nonparametric 
    % When to use Friedman (nonparametric):
    %   If differences are nonnormal or you prefer distribution‑free test: 
    %   run Friedman's test and follow with pairwise signrank tests with 
    %   multiple comparison correction.
    fprintf('\n --> Non-parametric. \n \n');
    p = friedman(mat,1,'off');
    fprintf('Friedman: \n   p-value = %1.4f \n',p);
    % If p_f < alpha: post-hoc pairwise with signrank + multiple comparison correction
    if p < alpha
        pairs = nchoosek(1:k, 2);
        rawP = zeros(k,1);
        for a = 1:k
            rawP(a) = signrank(mat(:,pairs(a,1)), mat(:,pairs(a,2)));
        end
        adjP = min(rawP * numel(rawP), 1); % Bonferroni
        fprintf('\nPaired comparisons with correction: \n');
        fprintf('\n \n   Raw p: %s\n   Adj p (Bonferroni): %s \n \n \n', ...
            mat2str(round(rawP,3)'), mat2str(round(adjP,3)'));
    end
end

%%
isParametric = ~all(normality) && sphericity > alpha;

switch isParametric
    case 1
        isParametric = 'true';
    case 0
        isParametric = 'false';
end