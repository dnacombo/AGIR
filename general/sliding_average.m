function y = sliding_average(x, dim, win, ic)
%  SLIDING_AVERAGE
%
%  Usage:
%    >> y = sliding_average(x, dim, win, ic);
%
%  Arguments:
%      x - input data
%    dim - dimension along which the sliding average is computed
%    win - average window (MUST BE AN ODD INTEGER)
%     ic - samples around which the sliding average is computed
%      y - averaged data

%  2007/12/14, Valentin Wyart (valentin.wyart@chups.jussieu.fr): v1.0
if not(exist('ic','var'))|| isempty(ic)
    ic = 1:size(x,dim);
end

dims = 1:ndims(x);

y = permute(x, [dim setdiff(dims, dim)]);
y_size = size(y);

y = filter(ones([1 win])/win, 1, y, [], 1);
y = y(min(floor(win/2)+ic, y_size(1)),:);
y = reshape(y, [length(ic) y_size(2:end)]);

y = ipermute(y, [dim setdiff(dims, dim)]);
