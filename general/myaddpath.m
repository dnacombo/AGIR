function myaddpath(varargin)

% myaddpath(varargin)
%
% same as addpath but converts relative paths to absolute paths.
%
% see also: addpath



i = 1;
while i <= numel(varargin)
    if strcmp(varargin{i}(1),'-')
        % option argument
        i = i+1;
        continue
    elseif not(isempty(regexp(varargin{i}(1:2),'[a-zA-Z]:', 'once')))
        % absolute path windows style
        i = i+1;
        continue
    elseif strcmp(varargin{i}(1),'/')
        % absolute path unix style
        i = i+1;
        continue
    elseif exist(varargin{i},'dir')
        % create absolute path
        varargin{i} = fullfile(pwd,varargin{i});
        i = i+1;
        continue
    end
end
addpath(varargin{:});