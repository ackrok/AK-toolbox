function [main] = shadedband(x, y)
% shadedband
%
% Description: This code creates a shaded vertical band from t = x1 to x2
%
% [main] = shadedband(x, y)
%
% Inputs: 'x' is 1x2 matrix with lower and upper limit along x-axis
%         'y' is same for y-axis
%
% Example usage:
%   figure;
%   shadedband([10, 11], ylim); 
% Will create shaded vertical shaded band from 10 to 11 spanning figure.
%
% Anya Krok, March 2026

xv = [x(1) x(2) x(2) x(1)]; 
yv = [y(1) y(1) y(2) y(2)];
h = patch(xv, yv, ...
    [0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.8, ...
    'DisplayName', '');
% uistack(h,'bottom'); % keep patch behind plotted data

end