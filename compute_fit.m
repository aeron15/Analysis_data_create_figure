function [fittedX, fittedY]=compute_fit(x,y)

%COMPUTE_FIT computes linear fit of x versus y

%% Removes NANs rom the data

matrix=[x' y'];
matrix=matrix(~any(isnan(matrix),2),:);
x=matrix(:,1);
y=matrix(:,2);

%%
coeffs = polyfit(x, y, 1);
% Get fitted values
fittedX = linspace(min(x)-1, max(x)+1, 200);
fittedY = polyval(coeffs, fittedX);
