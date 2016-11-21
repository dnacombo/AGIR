function [mid] = middle(x)

% [mid] = middle(x)
% returns the values that are between each value of x so that
% for all i in x, mid(i) = (x(i)+x(i+1))/2

mid = (x(1:end-1) + x(2:end)) /2;


