function [main] = shadedband(x, y, label)
% shadedband
% This code creates a shaded vertical band from t = x1 to x2.
%
% Syntax:
%   [main] = shadedband(x, y);
%   [main] = shadedband(x, y, label);
%
% Inputs: 'x' is 1x2 matrix with lower and upper limit along x-axis
%         'y' is same for y-axis
%         'label' (optional) char or string used as patch DisplayName
%                           (default is empty)
%
% Example usage:
%   figure;
%   shadedband([10, 11], ylim); 
% Will create shaded vertical shaded band from 10 to 11 spanning figure.
%
% Anya Krok, March 2026

if nargin < 3
    label = '';
end
if isstring(label); label = char(label); end

xv = [x(1) x(2) x(2) x(1)]; 
yv = [y(1) y(1) y(2) y(2)];
h = patch(xv, yv, ...
    [0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.8, ...
    'DisplayName', label);
% uistack(h,'bottom'); % keep patch behind plotted data

end