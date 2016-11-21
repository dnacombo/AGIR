function map = lines10(n)
%LINES10  Color map with the line colors.
%   LINES10(M) returns an M-by-3 matrix containing a "ColorOrder"
%   colormap. LINES10, by itself, is the same length as the current
%   colormap.
%   See also HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, 
%   COLORMAP, RGBPLOT.

if nargin<1, n = size(get(gcf,'Colormap'),1); end


c = [0,0,1
    0,0.5,0
    .8 0 0
    0,0.75,0.75
    0.75,0,0.75
    0.75,0.75,0
    0.3 .3 .3
    0,.7,0
    0.5,0.5,0
    0,0.5,0.5];


map = c(rem(0:n-1,size(c,1))+1,:);



