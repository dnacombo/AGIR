function fout = waitbarr(x,whichbar, varargin)

%WAITBARR Display wait bar * and estimated remaining time *.
%   H = WAITBAR(X,'message', property, value, property, value, ...)
%   creates and displays a waitbar of fractional length X.  The
%   handle to the waitbar figure is returned in H.
%   X should be between 0 and 1.  Optional arguments property and
%   value allow to set corresponding waitbar figure properties.
%   Property can also be an action keyword 'CreateCancelBtn', in
%   which case a cancel button will be added to the figure, and
%   the passed value string will be executed upon clicking on the
%   cancel button or the close figure button.
%
%   WAITBAR(X) will set the length of the bar in the most recently
%   created waitbar window to the fractional length X.
%
%   WAITBAR(X,H) will set the length of the bar in waitbar H
%   to the fractional length X.
%
%   WAITBAR(X,H,'message') will update the message text in
%   the waitbar figure, in addition to setting the fractional
%   length to X.
%
%   WAITBAR is typically used inside a FOR loop that performs a
%   lengthy computation.
%
%   Example:
%       h = waitbar(0,'Please wait...');
%       for i=1:1000,
%           % computation here %
%           waitbar(i/1000,h)
%       end
%
%   See also DIALOG, MSGBOX.

persistent t1 t2 x0 hs
if not(exist('x','var')) || isempty(x)
    x = 0;
end
msg = {};
if not(exist('whichbar', 'var'))
    whichbar = [];
elseif ischar(whichbar)
    msg = {whichbar};
    whichbar = [];
end
if ~isempty(hs)
    todel = [];
    for i = 1:numel(hs)
        if not(ishandle(hs(i)))
            todel(end+1) = i;
        end
    end
    hs(todel) = [];
    if not(isempty(whichbar))
        ibar = find(hs == whichbar);
    else
        ibar = numel(hs)+1;
    end
else
    ibar = 1;t1 = [];t2 = []; x0 = [];
end
est = Inf; 
if not(numel(t1) >= ibar)
    t1(ibar) = -5;% we'll run waitbarr 5 times before we start estimating time
elseif t1(ibar) < 0
    t1(ibar) = t1(ibar) +1;
elseif numel(x0) >= ibar && x0(ibar) > 0
    t2(ibar) = toc(uint64(t1(ibar)));
    est = 1 / (x - x0(ibar)) * t2(ibar) - t2(ibar);
else % when t1 reaches 0 we first come here, then 2 lines above
    t1(ibar) = tic; x0(ibar) = x;
end
if isinf(est)
    msg = {msg{:} [datestr(now)]};
else
    msg = {msg{:} [datestr(now)] [s2hms(round(est)) ' remaining']};
end
if isempty(whichbar)
    fout = waitbar(x,msg);
else
    fout = waitbar(x,whichbar,msg);
end
hs(ibar) = fout;
