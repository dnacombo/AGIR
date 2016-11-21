function barplot(data,varargin)
% barplot(data,commands)
% plot the average of data along dimension one and adds errorbars
% commands can be one of
% 'std' (default) error bars are standard deviation
% 'sem' to display standard error of the mean instead of stdev
% 'indiv' to overlay a scatterplot of individual datapoints of the
%           first dimension
% 'indivline' to overlay a scatterplot of individual datapoints of the
%           first dimension plus lines joining same indices along first dim
% 'indivlinelegend' also put a legend
% 'stars' to perform some basic ttest between columns of data and show
%           stars where paired differences are significant.
%
% Max 2012

cmds = varargin;
s = size(data);
if numel(s) > 2
    means = reshape(nanmean(data,1),s(2:end));
    stds = reshape(nanstd(data,[],1),s(2:end));
elseif numel(s) == 2
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
    if nBarsPerGroup == 1
        xl(2) = xl(2)-1;
    end
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
set(errorbarHandles,'tag','errorbars');
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
if numel(s) == 2
    data = reshape(data,[s(1) 1 s(2)]);
end
if not(isempty(regexpcell(cmds,'indiv')))
    ctrs = get(errorbarHandles,'Xdata');
    if iscell(ctrs)
        for i = 1:numel(ctrs)
            ct(:,i) = ctrs{i};
        end
    else
        ct = ctrs;
    end
    ctrs = permute(repmat(ct,[1,1,s(1)]),[3 1 2]);
    if not(isempty(regexpcell(cmds,'line')))
        %         for i_g = 1:nBarsPerGroup
        for i_d = 1:s(1)
            co = get(gca,'colororder');
            h(i_d,:) = plot(squeeze(ctrs(i_d,:,:))',squeeze(data(i_d,:,:))','o:','color',co(rem(i_d,size(co,1))+1,:));
            set(h(i_d,:),'tag',['indiv' num2str(i_d,'%02d')])
        end
            if not(isempty(regexpcell(cmds,'legend')))
                legend(h(:,1),'location','best')
            end
    else
        for i_bar = 1:nBarsPerGroup
            scatter(ctrs(:,i_bar),data(:,i_bar));
        end
    end
end
if not(isempty(regexpcell(cmds,'stars')))
    if numel(s) == 2
        yl = ylim;ystep = 0;
        for i_lim = 1:size(data,3)
            for i_lilim = i_lim+1:size(data,3)
                [h p ci stat] = ttest(data(:,:,i_lim),data(:,:,i_lilim));
                stars = pvalue2stars(p);
                if not(isnan(h))% && h
                    ystep = ystep + diff(yl)/20;
                    if all(get(hbar,'baseValue') > get(hbar,'yData'))
                        yl = yl(end:-1:1);
                    end
                    line([i_lim,i_lilim],[yl(2)-ystep yl(2)-ystep]);
                    text(mean([i_lim i_lilim]),yl(2)-ystep,stars,'horizontalalignment','center');
                end
            end
        end
    elseif numel(s) == 3
        yl = ylim;ystep = zeros(1,size(data,2));
        for i_lim = 1:size(data,3)
            for i_lilim = i_lim+1:size(data,3)
                [h p ci stat] = ttest(data(:,:,i_lim),data(:,:,i_lilim));
                stars = pvalue2stars(p);
                for i_h = 1:numel(h)
                    if not(isnan(h(i_h))) && h(i_h)
                        ystep(i_h) = ystep(i_h) + diff(yl)/20;
                        line([ctrs(i_lim,i_h,1),ctrs(i_lilim,i_h,1)],[yl(2)-ystep yl(2)-ystep]);
                        text(mean([ctrs(i_lim,i_h,1),ctrs(i_lilim,i_h,1)]),yl(2)-ystep,stars{i_h},'horizontalalignment','center');
                    end
                end
            end
        end
        
    end
end

if ~onhold
    hold off
end


function s = pvalue2stars(p)
% s = pvalue2stars(p)
%
% turn pvalue p (numeric) to cell array of strings s, the same size as p
% 
% p(i)<=0.05 : s{i} = '*'
% p(i)<0.01 : s{i} = '**'
% p(i)<0.001 : s{i} = '***'
% p(i)>0.5 : s{i} = 'NS'
% p(i)=NaN : s{i} = 'NaN'

% Proudly: Max August 2012

s = cell(size(p));
for i_p = 1:numel(p)
    
    if isnan(p(i_p))
        s{i_p} = 'NaN';
        continue
    end
    
    if p(i_p) > .05
        s{i_p} = 'NS';
    elseif p(i_p) < 0.001
        s{i_p} = '***';
    elseif p(i_p) < .01
        s{i_p} = '**';
    elseif p(i_p) <= .05
        s{i_p} = '*';
    end
end

