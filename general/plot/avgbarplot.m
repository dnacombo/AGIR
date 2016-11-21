function barplot(data,varargin)

% barplot(data,commands)
% plot the average of data along dimension one and adds errorbars (SEM)
% commands can be one of
% 'std' (default) error bars are standard deviation
% 'sem' to display standard error of the mean instead of stdev

cmds = varargin;
s = size(data);
if numel(s) > 2
    means = reshape(nanmean(data,1),s(2:end));
    stds = reshape(nanstd(data,[],1),s(2:end));
else
    means = nanmean(data,1);
    stds = nanstd(data,[],1);
end

if not(isempty(regexpcell(cmds,'sem')))
    stds = stds/sqrt(s(1));
end

[hbar] = bar(means);
if numel(size(means)) == 2
    xl = xlim;
    xlim([xl(1),ceil(xl(2)) - xl(1)])
end
onhold = ishold;
hold on
[errorbarHandles] = errorbar(means,stds,'LineStyle','none','color','black');
if ~onhold
    hold off
end

[nGroups, nBarsPerGroup] = size(means);
if nGroups == 1
    return
end
% The crazy stuff to find the bar centers. Some of it is how bar.m does it.
% This comes from barerrorbar
barHandlesStruct = get(hbar);
barWidth = barHandlesStruct(1).BarWidth; 

for iBpg = 1:nBarsPerGroup
   groupMembership(:,iBpg) = barHandlesStruct(iBpg).XData;
end
groupWidth = min(0.8, nBarsPerGroup/(nBarsPerGroup+1.5));
groupLocs = groupMembership(:,1);
distanceToNearestGroup = zeros(1,nGroups);
for iGroup = 1:nGroups
    distanceToNearestGroup(iGroup) = ...
        min(abs(groupLocs(iGroup)-...
        groupLocs(groupLocs~=groupLocs(iGroup))));
end
groupWidth = groupWidth*min(distanceToNearestGroup);
barGap = (nGroups - 1) * groupWidth / (nGroups-1) / nBarsPerGroup;
almostCenters  = (0:nBarsPerGroup-1)'.*barGap - 0.5*barGap*nBarsPerGroup;
relativeCenters = almostCenters + mean([(1-barWidth)/2.*barGap; (1+barWidth)/2.*barGap]);

centers = repmat(relativeCenters',nGroups,1) + groupMembership;
% Change the XData of the errorbars to be at our bar centers.
for iBpg = 1:nBarsPerGroup
    set(errorbarHandles(iBpg),'XData',centers(:,iBpg));
end

