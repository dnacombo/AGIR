function Xo = limo_orth(X)

% Xo = limo_orth(X)
% orthogonalization of each column of X with respect to the previous ones in sequence.
% this function just uses spm_orth. Displays the correlation matrix of X
% before and after orthogonalization, just for information.
%
% see also spm_orth.m



Xo = spm_orth(X);

disp('correlation between regressors after orthogonalization.')
[rho pval] = corr(Xo)





