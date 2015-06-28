function [fittedX, fittedY]=compute_fit(x,y)

%COMPUTE_FIT computes fit of x versus y

coeffs = polyfit(x, y, 1);
% Get fitted values
fittedX = linspace(min(x)-1, max(x)+1, 200);
fittedY = polyval(coeffs, fittedX);
