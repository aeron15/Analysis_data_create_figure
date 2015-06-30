function [x, y]=remove_nan_rows(x,y)

%REMOVE_NAN_ROWS removes NAN rows for correlation computations

%% Removes NANs rom the data

matrix=[x' y'];
matrix=matrix(~any(isnan(matrix),2),:);
x=matrix(:,1);
y=matrix(:,2);

end