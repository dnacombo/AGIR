function avgplot(data,varargin)

% barplot(data,commands)
% plot the average of data along dimension one with errorbars (SEM)
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

onhold = ishold;
hold on
[errorbarHandles] = errorbar(means,stds);
if ~onhold
    hold off
end


