function [varargout] = cell2coma(c)

% [comaseparatedlist] = pass(c)
% return the values of each cell of c as a coma separated list.
% 
% NB: this is equivalent to [comaseparatedlist] = c{:};

for i = 1:nargout
    varargout{i} = c{i};
end
