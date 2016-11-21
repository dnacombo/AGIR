function errorbarHandles = stackedplot(data,varargin)
% stackedplot(data,commands)
% plot the average of data along dimension one and adds errorbars...
% just look what it does... data should be 2 or 3 D.
% commands can be one of
% 'std' (default) error bars are standard deviation
% 'sem' to display standard error of the mean instead of stdev
% 'indiv' to overlay a scatterplot of individual datapoints of the
%           first dimension
% 'indivline' to overlay a scatterplot of individual datapoints of the
%           first dimension plus lines joining same indices along first dim
% 'indivlinelegend' also put a legend
%
% Max 2012

cmds = varargin;
s = size(data);
if numel(s) > 2
    means = reshape(nanmean(data,1),s(2:end));
    stds = reshape(nanstd(data,[],1),s(2:end));
else
    means = nanmean(data,1)';
    stds = nanstd(data,[],1)';
end

if not(isempty(regexpcell(cmds,'sem')))
    stds = stds/sqrt(s(1));
end

[ nBarsPerGroup, nGroups] = size(means);

% here's a trick for when we have only one value along dimension 1.
% Otherwise bar would display bars in one color.
if nBarsPerGroup == 1 % if there's only one
    means = cat(1,means,NaN(size(means)));
    [hbar] = bar(means);
    means(2,:) = [];
else
    [hbar] = bar(means);
end
if numel(size(means)) == 2
    xl = xlim;
%     if nBarsPerGroup == 1
%         xl(2) = xl(2)-1;
%     end
    xlim([xl(1),ceil(xl(2)) - xl(1)])
end

set(hbar,'tag','bars');
onhold = ishold;
hold on
if nBarsPerGroup == 1
    means = cat(1,means,NaN(size(means)));
    stds =  cat(1,stds,NaN(size(stds)));
end
[errorbarHandles] = errorbar(means,stds,'LineStyle','none','color','black');
%set(errorbarHandles,'tag','errorbars');
if nBarsPerGroup == 1
    means(2,:) = [];
    stds(2,:) = [];
end


barHandlesStruct = get(hbar);
% Change the XData of the errorbars to be at our bar centers.
if nGroups > 1
    for i_G = 1:nGroups
        vtc = get(get(hbar(i_G),'Children'),'Vertices');
        goodidx = logical([repmat([0 0 1 1 0]',nBarsPerGroup,1); 0]);
        ctrs = mean(reshape(vtc(goodidx,1),2,nBarsPerGroup));
        if nBarsPerGroup == 1
            ctrs = [ctrs NaN];
        end
        set(errorbarHandles(i_G),'XData', ctrs)
    end
end
delete(hbar);
xctrs = get(errorbarHandles,'Xdata');
yctrs = get(errorbarHandles,'Ydata');
uctrs = get(errorbarHandles,'Udata');
lctrs = get(errorbarHandles,'Ldata');
if iscell(xctrs)
    for i = 1:numel(xctrs)
        xct(:,i) = xctrs{i};
    end
else
    xct = xctrs;
end
if iscell(yctrs)
    for i = 1:numel(yctrs)
        yct(:,i) = yctrs{i};
    end
else
    yct = yctrs;
end
if iscell(uctrs)
    for i = 1:numel(uctrs)
        uct(:,i) = uctrs{i};
    end
else
    uct = uctrs;
end
if iscell(lctrs)
    for i = 1:numel(lctrs)
        lct(:,i) = lctrs{i};
    end
else
    lct = lctrs;
end
delete(errorbarHandles);
errorbarHandles = errorbar(xct',yct',lct',uct','color','k');
set(errorbarHandles,'tag','errorbars');
if not(isempty(regexpcell(cmds,'indiv')))
    ctrs = permute(repmat(xct,[1,1,s(1)]),[3 1 2]);
    if not(isempty(regexpcell(cmds,'line')))
        if numel(s) == 2
            data = reshape(data,[s(1) 1 s(2)]);
        end
        %         for i_g = 1:nBarsPerGroup
        for i_d = 1:s(1)
            co = get(gca,'colororder');
            h(i_d,:) = plot(squeeze(ctrs(i_d,:,:)),squeeze(data(i_d,:,:)),'o:','color',co(rem(i_d,size(co,1))+1,:));
            set(h(i_d,:),'tag',['indiv' num2str(i_d,'%02d')])
        end
        %         if i_g == 1
        if not(isempty(regexpcell(cmds,'legend')))
            legend(h(:,1),'location','best')
        end
        %         end
        %         end
    else
        for i_bar = 1:nBarsPerGroup
            scatter(ctrs(:,i_bar),data(:,i_bar));
        end
    end
end

if ~onhold
    hold off
end







