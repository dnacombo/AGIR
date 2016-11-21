function [p r c] = subpanel(r,c)

% create a panel with a layout specified in 1 of 2 ways:
% 
% [p r c] = subpanel(n)
% creates an approximately square panel
% later select it with p(row,col).select();
%
% [p r c] = subpanel(r,c)
% creates a panel with r rows and c columns

    
if not(exist('c','var'))
    [r c] = num2rowcol(r);
end

clf;
p = panel();
p.pack(r,c);


