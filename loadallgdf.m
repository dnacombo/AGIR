function [fs,x,event,C] = loadallgdf(datadir,participant)


fs = dir(fullfile(datadir,participant,'*.gdf'));% list all gdf files
fs = {fs.name};
fs(~emptycells(strfind(fs,'S7'))) = [];% remove S7

%%%%%%%
%%%%%%%%
load(fullfile(datadir,participant,['Task_' regexprep(participant,'(\d{1,2})','_$1')]));
list_trig = trig_theoretical(C);
%%%%%%%%
%%%%%%%%
data = {};
for i_f = 1:numel(fs)
    
    
    %%%%%%%
    datafile = fullfile(datadir,participant,fs{i_f});
    %%%%%%%%
    
    cfg = [];
    cfg.dataset = datafile;
    cfg.continuous = 'yes';
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [49 51];
    data{i_f} = ft_preprocessing(cfg);
    data{i_f}.event = ft_read_event(cfg.dataset);
end
%%

alldata = ft_appenddata([],data{:});
alldata.event = [];
for i_f = 1:numel(data)
    evt = data{i_f}.event;
    [evt.sample] = rep2struct([evt.sample] + sum(cellfun(@(x)size(x.trial{1},2),data(1:i_f-1))));
    alldata.event = [alldata.event evt];
end
[alldata.event.value] = rep2struct([alldata.event.value] - 65280);
if not(ismember(250,[alldata.event.value]))
    [alldata.event.value] = rep2struct([alldata.event.value] - 1);
end
if not(ismember(250,[alldata.event.value]))
    figure;scatter(1:numel(alldata.event),[alldata.event.value])
    error('deal with triggers codes')
end

event = alldata.event;
% figure;hist([event.value],1:255)
% xtick(unique([event.value]))
% xticklabel(unique([event.value]))

tb = struct2table(event);
tb = tb([tb.value == 30] | tb.value == 68 | tb.value == 240,:);
list_trig = list_trig([list_trig == 30] | list_trig == 68 | list_trig == 240 | list_trig == 240);
if not(isequal(tb.value',list_trig))
    disp('Warning: some triggers missing?!? Attempting a repair...')
    
    
end

fs = alldata.fsample;

x = [alldata.trial{:}];
