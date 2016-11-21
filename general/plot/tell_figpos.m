function tell_figpos(where)

% tell_figpos(where)
% display figure position stored in global variable figpos to command window
% if where is 'all', display all stored figure positions
%
% see also copy_figpos and paste_figpos

global figpos
if nargin == 1
    if strcmp(where,'all')
        disp(getpref('figpos'))
        return
    end
    try
        figpos = getpref('figpos',where);
    end
end
if not(isempty(figpos))
    disp(figpos);
end