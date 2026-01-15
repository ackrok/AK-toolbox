function tf = isVecInteger(vec, tol)
% Returns true if v (double vector/array) has all finite elements integer
% within tol. NaNs are ignored/treated as acceptable. Inf values fail.
%
% tf = isVecInteger(vec)
% tf = isVecInteger(vec, tol)
%
% INPUTS
% 'vec' - vector of values
% 'tol' - tolerance (optional), default is 0
%       Use tol = 0 for exact checks.
%       Use a small positive tol to tolerate floating‑point roundoff.
% OUTPUTS
% 'tf' - true/false return if vector has finite integer elements.
%
% Anya Krok, January 2026

if nargin < 2, tol = 0; end
validateattributes(vec, {'double','single','numeric'}, {}, mfilename, 'v', 1);
tf = isvector(vec) && all( ~isfinite(vec) | (isnan(vec) | abs(vec - round(vec)) <= tol) );
% Above: isfinite(v) false for NaN/Inf. We want to
%  - ignore NaN (treat as OK)
%  - reject Inf
% So keep final check clearer below:
finiteMask = isfinite(vec);
if any(isinf(vec))
    tf = false;
    return
end
tf = all(abs(vec(finiteMask) - round(vec(finiteMask))) <= tol);
end