% uigetfile2() - same as uigetfile but remember folder location.
%
% Usage: >> uigetfile2(...)
%
% Inputs: Same as uigetfile
%
% Author: Arnaud Delorme & Hilit Serby, Scott Makeig, SCCN, UCSD, 2004
%         Thanks to input from Bas Kortmann
% Modified to use user specific matlab preference utility, Maximilien Chaumon, 2011 
%
%
% Copyright (C) Arnaud Delorme & Hilit Serby, Scott Makeig, SCCN, UCSD, 2004

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function varargout = uigetfile2(varargin);
mlock
persistent tmp_fld

if nargin < 1
    help uigetfile2;
    return;
end;

% remember old folder
olddir = pwd;
try
    eeglab_options;
    if option_rememberfolder
        if isempty(tmp_fld)
            tmp_fld = getpref('eeglab','lastdir',olddir);
        end
        cd(tmp_fld);
    end
end

%% Show the open dialog and save the latest directory to the file
% ---------------------------------------------------------------
[varargout{1} varargout{2} varargout{3}] = uigetfile(varargin{:});
try
    if option_rememberfolder
        if ischar(varargout{1})
            tmp_fld = varargout{2};
            setpref('eeglab','lastdir',tmp_fld);
        end
    end
end
cd(olddir)

