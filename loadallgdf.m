function [fs,x,event,C] = loadallgdf(datadir,participant)


load(fullfile(datadir,participant,['Task_' regexprep(participant,'(\d{1,2})','_$1')]));

lastevent = 220;

fs = dir(fullfile(datadir,participant,'*.gdf'));% list all gdf files
fs = {fs.name};
switch participant
    case 'participant19'
        fs(regexpcell(fs,'S1|S2')) = [];% remove S1 S2 which were defect
        % remove blocks 15 to 17
        C(15:17) = [];
    case 'participant22'
        fs(regexpcell(fs,'x')) = [];% remove two _x_ files which were defect
        % remove blocks 55 to last
        C(54:end) = [];
    otherwise
        fs(~emptycells(strfind(fs,'S7'))) = [];% remove S7
end

%%%%%%%
%%%%%%%%
expected_triggs = trig_theoretical(C);
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
if strcmp(participant,'participant22')
    event = event(1:10045);
    event(end).value = 220;
end
% figure;hist([event.value],1:255)
% xtick(unique([event.value]))
% xticklabel(unique([event.value]))

allevents = struct2table(event);
received_triggs = allevents([allevents.value == 30] | allevents.value == 68 | allevents.value == 240,:);
expected_triggs = expected_triggs([expected_triggs == 30] | expected_triggs == 68 | expected_triggs == 240);
if size(received_triggs,1) > numel(expected_triggs)
    last = find(allevents.value==lastevent,1,'last');
    received_triggs = received_triggs(1:numel(expected_triggs),:);
    allevents = allevents(1:last,:);
end
    
if not(isequal(received_triggs.value',expected_triggs))
    col_sync=2;
    i_received = 1;
    i_block = 0;
    currsync = 68; % by default premier run est sync
    syncasync = [68,30];% corresponding in C: 1=sync, 2=async
    while i_received <= numel(received_triggs.value')
        if ismember(received_triggs.value(i_received),[30 68])
            currsync = received_triggs.value(i_received);
        end
        if ismember(received_triggs.value(i_received),[240])
            i_block = i_block+1;
            C{i_block}(:,9) = 1;% block is ok a priori. If it's bad, we'll change this below
        end
        if received_triggs.value(i_received) ~= expected_triggs(i_received)
            i_val_allevents = allevents.sample == received_triggs.sample(i_received);
            switch expected_triggs(i_received)
                case {30, 68}
                    % list_trig_synty_def = 30/68 (définition du bloc comme synch / asynch)
                    % Utilisée pour définir les blocs + pour exclure ceux où le manque du trigger fait que le script run un bloc inattendu.
                    % - Récupérer définition attendue de [ce bloc],
                    % comparer à ce qu'on à sur la base des triggs reçus
                    % jusqu'ici
                    shouldbesync = syncasync(C{i_block}(1,col_sync));
                    if currsync ~= shouldbesync
                        C{i_block}(:,9) = 0;
                    end
                    toadd = allevents(i_val_allevents,:);
                    toadd.type = {'AJOUTE'};
                    toadd.value = currsync;
                    allevents(end+1,:) = toadd;
                    %                 case 240
                    % this actually should never happen, we're missing only
                    % 30s and 68s...
                    %                     % Utilisée seulement pour définir les blocs.
                    %                     % Repair: soit ajouter un trigger artificiel entre les deux triggers voisins (en amont 250 ou 80: fin de l’échelle ou début de session, et en aval 112/150: premier flip du bloc); soit juste remplacer par le premier flip du bloc (ou premier flip du bloc .sample-1 pour ne pas que ce soit exactement les mêmes..?)...
                    %                     all_future_events =  allevents(i_val_allevents_last:end,:);
                    %                     all_future_112_150 = all_future_events(ismember(all_future_events.value,[112 150]),:);
                    %                     first_flip_next_block = all_future_112_150(1,:);
                    %                     first_flip_next_block.type = {'AJOUTE'};
%                     first_flip_next_block.value = expected_triggs(i_val);
%                     allevents(end+1,:) = first_flip_next_block;
                otherwise
                    error('unexpected')
            end
            expected_triggs(i_received) = [];
        end
        i_received = i_received +1;
    end
end
if not(isequal(expected_triggs',received_triggs.value))
    error('messed up')
end
event = table2struct(sortrows(allevents,'sample'));
fs = alldata.fsample;

x = [alldata.trial{:}];
