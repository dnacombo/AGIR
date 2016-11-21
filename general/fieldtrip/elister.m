function [varargout] = elister(re,evts,varargin)

% [events] = elister(re,evts)
% [fstruct] = elister(re,evts,'key', val, ...)
% 
% This function lists all events matching re in the event structure evts.
% re is a regular expression with named tokens
%
% optional input key val pairs : 
%           'exclude', str  : file names matching regexp str will be
%                            excluded from the list 
%           'dir',str       : root directory where to search for the files 
% outputs:
%           cellstrevts: cell of events
%           factnames: cell: all the factor names of the named tokens of re
%           factlevelnames cell: the level names of each factor
%           factidx:   cell: for each factor, the level to which each file belongs
% output method 2:
%           f: 	structure with fields name (name of the file), and each of
%               the factor names with the corresponding level to which the file belongs.
if numel(varargin) == 1 && isstruct(varargin{1})
    varargin = struct2vararg(varargin{1});
end
g = finputcheck( varargin, ...
    {
    'exclude'   { 'string';'cell' }    []   ''
    'dir' { 'string'}    []   cd
    });
rootdir = g.dir;

allfields = fieldnames(evts);
evtsfields = regexp(re,'<(.*?)>','tokens');
for i = 1:numel(evtsfields)
    evtsfields{i} = evtsfields{i}{1};
end
evtsfields = allfields(regexpcell(allfields,evtsfields,'exact'));

cellstrevts = cell(numel(evts),1);
for i_ev = 1:numel(evts)
    str = ['event:' num2str(i_ev)];
    for i_f = 1:numel(evtsfields)
        if isnumeric(evts(i_ev).(evtsfields{i_f}))
            str = [str ';' evtsfields{i_f} ':' num2str(evts(i_ev).(evtsfields{i_f})) ];
        elseif ischar(evts(i_ev).(evtsfields{i_f}))
            str = [str ';' evtsfields{i_f} ':' evts(i_ev).(evtsfields{i_f}) ];
        elseif iscellstr(evts(i_ev).(evtsfields{i_f}))
            str = [str ';' evtsfields{i_f} ':' evts(i_ev).(evtsfields{i_f}){1} ];
        end
        if i_f == numel(evtsfields)
            str = [str ';'];
        end
    end
    cellstrevts{i_ev} = str;
end
names = regexp(cellstrevts,re,'names');
if isempty(names) || numel(fieldnames(names{1})) == 0
    disp('Warning: Could not find tokens in regular expression')
    factnames = {};
    f = struct('name',cellstrevts);
else
    todel = [];
    for i = 1:numel(names)
        if isempty(names{i})
            todel(end+1) = i;% removing events that did not match the whole pattern
        else
            f(i) = names{i};
        end
    end
    if isempty(setxor(todel, 1:numel(names)))
        error('Could not find anything. Check your search pattern.')
    end
    f(todel(todel<=numel(f))) = [];
    cellstrevts(todel) = [];
    factnames = fieldnames(f);
    [f.name] = cellstrevts{:};
end
todel = [];
for i = 1:numel(factnames)
    factlevelnames{i} = unique({f.(factnames{i})});
    if numel(factlevelnames{i}) == 1
        todel(end+1) = i; 
        % remove uninformative fields, that have only one value and in every trial
    end
    finder = regexptranslate('escape',{f.(factnames{i})});
    factidx{i} = regexpcell(factlevelnames{i},finder,'exact');
end
f = rmfield(f,factnames(todel));% remove here
for i = 1:numel(cellstrevts)
    for j = 1:numel(factnames)
        if not(any(j == todel))% and don't add here
            f(i).([factnames{j} 'idx']) = factidx{j}(i);
        end
    end
end

if nargout > 1
    varargout{1} = cellstrevts;
    varargout{2} = factnames;
    varargout{3} = factidx;
    varargout{4} = factlevelnames;
else
    varargout{1} = f;
end