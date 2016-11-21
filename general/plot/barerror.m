function barerror(data, varargin)
% function barerror(data, varargin)
% Plots bars for each row with error bars computed as range of data
%
% USAGE:
% Pass in matrix of data arranged with repeated measurements in columns,
% with rows containing individual types
%           -OR-
% A cell array of these data matrices, for plotting multiple adjacent bars with
% unique colors.
%
% The rest of the stuff it expects goes as follows:
% 1) method: either 1 for SEM, 2 for 95% CI, 3 for min max errorbars
% 2) barlabels:  cell array of strings to label the ROWS of data
% 3) grouplabels: cell array of strings to label the matrix/matrices of data 
% 4) colors: the sequence of color characters for labeling bars
%
% Jonathan Scholz,  GaTech, October 15 2009

% Example:
% data = {diag(sin(1:5)) * rand(5,5),diag(cos(1:5)) * rand(5,5)};
% barerror(data,{'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'},{'Sin','Cos'});

if iscell(data)
    ntypes = size(data,2);
    nbars = size(data{1},1);
else
    ntypes = 1;
    nbars = size(data,1);
end

if nargin >=2 && not(isempty(varargin{1}))
    method = varargin{1};
else
    method = 1; % default to SEM error bars
end

if nargin >= 3 && not(isempty(varargin{2}))
    barlabels = varargin{2};
else
    for i=1:nbars
        barlabels{i}=sprintf('%d',i);
    end
end

if nargin >= 4 && not(isempty(varargin{3}))
    grouplabels = varargin{3};
else
    for i=1:nbars
        grouplabels{i}=sprintf('Class %d',i-1);
    end
end

if nargin >= 5 && not(isempty(varargin{4}))
    colors = varargin{4};
else
    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'];
end

clf;
hold on;
barhandles=[];
for i=1:ntypes
    if iscell(data)
        d=data{i};
    else
        d=data;
    end

    % data stuff
    means = nanmean(d,2)';
    mins = min(d,[],2)';
    maxs = max(d,[],2)';
    stdev = nanstd(d,[],2)';
    sem = stdev/sqrt(size(d,2));
    if method == 1
        L = sem;
        U = sem;
    elseif method == 2
        Z = 1.96; % for 0.95 CI
        L = - Z*sem';
        U = Z*sem';
    elseif method == 3
        L = means-mins;
        U = means+maxs;
    end

    % plot control stuff
    xrange = 1:nbars;
    offsets = fix(-ntypes/2):fix(ntypes/2);
    smooshFactor = 0.8;
    barOrigins = xrange+smooshFactor*(offsets(i)/ntypes);
    barwidth = smooshFactor/ntypes;

    b = bar(barOrigins, means, barwidth, colors(mod(i,length(colors))));
    barhandles=[barhandles b]; % is there a better way to fix this???
    errorbar(barOrigins,means,L,U,'LineStyle','none','color','black');
end

set(gca,'XTickLabel',barlabels,'xtick',xrange);
% title_str = 'title';
% xlabel_str = 'x axis';
% ylabel_str = 'y axis';

% set(gca,'FontSize',16)
% 
% title(title_str);
% h = get(gca, 'title');
% set(h,'FontSize',16);
% 
% xlabel(xlabel_str);
% h = get(gca, 'xlabel');
% set(h,'FontSize',24);
% 
% ylabel(ylabel_str);
% h = get(gca, 'ylabel');
% set(h,'FontSize',24);

if ntypes > 1
    legend(barhandles, grouplabels);
end
