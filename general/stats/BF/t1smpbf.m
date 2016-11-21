function bf10 = t1smpbf(t,n,r)
%
% bf10 = t1smpbf(t,n,[r=0.707])
%
% Calculates JZS Bayes Factor for a one-sample t-test given t and sample size n.
% The optional input r is the scale factor which defaults to 0.707.
% This quantifies the evidence in favour of the alternative hypothesis. 
% See Rouder et al, 2009, Psychon Bull Rev for details.
%

% Default scale factor
if nargin < 3
    r = 0.707;
end

% Function to be integrated
F = @(g,t,n,r) (1+n.*g.*r.^2).^(-1./2) .* (1 + t.^2./((1+n.*g.*r.^2).*(n-1))).^(-n./2) .* (2.*pi).^(-1./2) .* g.^(-3./2) .* exp(-1./(2.*g));

% Bayes factor calculation
bf01 = NaN(size(t));
% note this is faster than using arrayfun
if numel(n) == 1
    n = repmat(n,size(t));
end
for i = 1:numel(t)
    bf01(i) = (1 + t(i)^2/(n(i)-1))^(-n(i)/2) / integral(@(g) F(g,t(i),n(i),r),0,Inf);
end

% Invert Bayes Factor
bf10 = 1 ./ bf01;