function [selx,sely,indexnr]=lasso(x,y,varargin)

% lasso -  enables the selection/encircling of (clusters of) events in a scatter plot by hand
%          using the mouse
%
% Input:    x,y                 - a set of points in 2 column vectors.
% Output:   selx,sely,indexnr   - a set of selected points in 3 column vectors
%
% Note:   After the scatter plot is given, selection by mouse is started after any key press.
%         This is done to be able to ZOOM or CHANGE AXES etc. in the representation before selection
%         by mouse.
%         Encircling is done by pressing subsequently the LEFT button mouse at the requested positions
%         in a scatter plot.
%         Closing the loop is done by a RIGHT button press.
%
% T.Rutten V2.0/9/2003

opts = vararg2struct(varargin);

h = figure;
plot(x,y,'k.')

if isfield(opts,'title')
    title(opts.title);
end
if isfield(opts,'axis')
    axis(opts.axis);
end

las_x=[];
las_y=[];

c=1;

key=0;

while c==1
    [a,b,c]=ginput(1);
    las_x=[las_x;a];las_y=[las_y;b];
    line(las_x,las_y)
    drawnow
end;

las_x(length(las_x)+1)=las_x(1);
las_y(length(las_y)+1)=las_y(1);

line(las_x,las_y)
pause(.2)

in=inpolygon(x,y,las_x,las_y);

ev_in=find(in>0);

selx=x(ev_in);
sely=y(ev_in);

figure(h);clf
plot(x,y,'b.',selx,sely,'g.');
legend(num2str([length(x)-length(selx);length(selx)]));
if isfield(opts,'axis')
    axis(opts.axis);
end

indexnr=ev_in;