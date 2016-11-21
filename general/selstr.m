function idx = selstr(c,pat,inv)

% idx = selstr(c,pat,inv)
% return indices of cells in c that match pattern pat (regexp). 
% invert if inv is true.
% pat can be char or cellstr
error('use regexpcell')

error(nargchk(2,3,nargin))
if not(iscellstr(c))
    error('input c must be a cell array of strings');
end

idx = [];
if ischar(pat)
    pat = cellstr(pat);
end

for i_pat = 1:length(pat)
    trouv = regexp(c,pat{i_pat});
    for i = 1:numel(trouv)
        if not(isempty(trouv{i}))
            idx(end+1) = i;
        end
    end
end
if exist('inv','var') && inv
    others = 1:numel(trouv);
    others(idx) = [];
    idx = others;
end