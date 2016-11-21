function r = randtruncexp(mu,trunc,siz)
%randtruncexp Random arrays from truncated exponential distribution.
%   R = randtruncexp(mu,trunc) returns an array of random numbers chosen from the
%   exponential distribution with mean parameter mu, Truncated at trunc.
%   The size of R is the size of MU.
%
%   R = randtruncexp(mu,trunc,siz) returns an array of size siz.
%
%   See also EXPCDF, EXPFIT, EXPINV, EXPLIKE, EXPPDF, EXPSTAT, RANDOM.

% Maximilien Chaumon, based on 
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/81616

if nargin < 2
    error('Not enough input arguments')
end
if nargin < 3
    siz = size(mu);
end

r = -log(1-rand([siz,1]).*(1-exp(-trunc./mu))).*mu;

