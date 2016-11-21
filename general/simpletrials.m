function T = simpletrials(conds,names,varargin)


for i_tri = 1:numel(conds{1})
    for i_cond = 1:numel(conds)
        T(i_tri).(names{i_cond}) = conds{i_cond}(i_tri);
    end
end