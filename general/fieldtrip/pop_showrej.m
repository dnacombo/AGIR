

function [chanlist strchannames cellchannames] = pop_showrej(EEG)

if nargin == 0
    EEG = evalin('base','EEG');
end

h0 = figure(23537);clf;
set(gcf,'numbertitle','off','name','Channel selection');
topoplot([],EEG.chanlocs,'electrodes','on', 'style','blank');

ax = get(gcf,'children');
chi = get(ax,'children');
chi(end-4:end) = [];
hold on
X = get(chi,'XData');
X = [X{:}];
Y = get(chi,'YData');
Y = [Y{:}];
delete(chi)
chi = get(ax,'children');
delete(chi(end-4));

set(gca,'position',[0 0 1 1]);% take the whole figure space for the axis
set(gca,'units','pixels')
pp = get(gca,'position');% size of the axis in pixels
set(gca,'units','normalized');

% now I'll try to make only one surface, because 64 is too slow.
sp = 10;% number of pixels per ERPimg
ss = sp./pp(3:4);% representing ss units in head space of X and Y
XX = linspace(min(X)-ss(1),max(X)+ss(1),pp(3)); % pixel space
YY = linspace(min(Y)-ss(2),max(Y)+ss(2),pp(4));

ZZ = zeros(numel(XX),numel(YY));
CC = NaN*zeros(numel(XX),numel(YY));
% now we would like to fill ZZ and CC appropriately
% ZZ must be NaNs, except where we have electrodes, where it'll be 1
% CC must be NaNs, except where we have electrodes, where it'll be EEG.data
for i = 1:numel(X)
    % x y are the indices of XX and YY where the electrodes should be
    % drawn.
    x = find(abs(XX - X(i)) == min(abs(XX - X(i))));
    y = find(abs(YY - Y(i)) == min(abs(YY - Y(i))));
    
    ZZ(x-sp:x+sp,y-sp:y+sp) = 1;% ones where we want to draw the elecs.
end
%     xx = X(i)-surfsizes:XXsr:X(i)+surfsizes;%linspace(X(i)-surfsizes,X(i)+surfsizes,EEG.pnts);
%     yy = Y(i)-surfsizes:YYsr:Y(i)+surfsizes;%linspace(Y(i)-surfsizes,Y(i)+surfsizes,EEG.trials);
%     zz = zeros(numel(xx),numel(yy))';
%     cc = double(squeeze(EEG.data(i,:,:))');
    h = surf(XX,YY,ZZ');
    shading flat;



