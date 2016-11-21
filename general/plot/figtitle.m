function figtitle(str,varargin)

% h = figtitle(str,varargin)
% create a title centered to the top of the current figure
% with the string (or cell array of strings) str
% optional argument pairs to uicontrol('style','text', ...)
% are passed via varargin

h = gcf;

%%%% delete previous figtitle if any
hh = findobj(h,'tag', 'FigTitle');
if ~isempty(hh)
    delete(hh);
end
%%%% keep consistent toolbar appearance
hastoolbar = get(h,'toolbar');
% if strcmp(hastoolbar,'auto')
%     hastoolbar = findobj(h,'type', 'uicontrol');
%     if not(isempty(hastoolbar))
%         hastoolbar = 'none';
%     else
%         hastoolbar = 'figure';
%     end
% end
hasmenubar = get(h,'menubar');
%%%%

%%%% create new figtitle
hh = uicontrol('style','text','backgroundcolor',get(h,'color'),...
    'units','normalized','position',[.5 1 .0001 .0001],'string',str,...
    'horizontalalignment','center','tag','FigTitle',varargin{:});
%%%% set proper position according to fontsize and number of lines 
set(hh,'units','characters')
p = get(hh,'position');
ff = get(hh,'fontsize')./12;
switch class(str)
    case 'char'
        w = size(str,2);
        set(hh,'position',p + ff*[-w -size(str,1)-1 2*w size(str,1)+1])
    case 'cell'
        w = max(cellfun(@numel,str));
        set(hh,'position',p + ff*[-w -numel(str)-1 2*w numel(str)+1])
end
set(hh,'units','normalize')
%%%%
set(gcf,'toolBar',hastoolbar)
set(gcf,'menubar',hasmenubar)
uistack(hh,'bottom');

