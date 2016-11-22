function [fs,x,event,C] = loadallgdf(datadir,participant)


fs = dir(fullfile(datadir,participant,'*.gdf'));% list all gdf files
fs = {fs.name};
fs(~emptycells(strfind(fs,'S7'))) = [];% remove S7

%%%%%%%
%%%%%%%%
load(fullfile(datadir,participant,['Task_' regexprep(participant,'(\d{1,2})','_$1')]));
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
% figure;hist([event.value],1:255)
% xtick(unique([event.value]))
% xticklabel(unique([event.value]))

allevents = struct2table(event);
received_triggs = allevents([allevents.value == 30] | allevents.value == 68 | allevents.value == 240,:);
expected_triggs = expected_triggs([expected_triggs == 30] | expected_triggs == 68 | expected_triggs == 240);
if not(isequal(received_triggs.value',expected_triggs))
    allC = vertcat(C{:});
    col_nb_total_trial=1;
    col_sync=2;
    col_percentage = 3;
    col_this_trial_type = 4;
    col_what_proportion_tells_type = 5;
    
    for i_val = 1:numel(received_triggs.value')
        % what is current block?
        if received_triggs.value(i_val) ~= expected_triggs(i_val)
            % previous trigger is the last one that was correct (i.e. for
            % which we have a sample value)
            last_received_trigger_sample = received_triggs.sample(i_val-1);
            % position in list of all events (find position of previous trigger):
            i_val_allevents_last = dsearchn(allevents.sample,last_received_trigger_sample);
            switch expected_triggs(i_val)
                case {30, 68}
                    % list_trig_synty_def = 30/68 (définition du bloc comme synch / asynch)
                    % Utilisée pour définir les blocs + pour exclure ceux où le manque du trigger fait que le script run un bloc inattendu.
                    % Repair:
                    % - Récupérer définition de [ce bloc] et [ce bloc -1] dans mat C
                    issync = allC(i_val,col_sync);
                    prev_issync = allC(i_val-1,col_sync);
                    % - si def_ce_bloc == bloc-1 :
                    if issync == prev_issync
                        % Soit rajouter un trigger de la définition correspondante juste avant le prochain trig de begin scale (180) soit se servir directement du trigger de begin scale.
                        allevents(i_val_allevents_last-10:i_val_allevents_last+10,:)% show surroundings
                        disp('There should be a 180 around here...')
                        keyboard
                    else
                        % - si ce bloc =! bloc d’avant : exclure le bloc (p.e. en enlevant tous les triggers jusqu’à la définition du bloc suivant ? Ou définir une sorte de red flag.)
                        all_future_events =  allevents(i_val_allevents_last:end,:);
                        all_future_block_begin = all_future_events(ismember(all_future_events.value,[30 68]),:);
                        next_block_begin = all_future_block_begin(1,:);
                        all_previous_events =  allevents(i_val_allevents_last:end,:);
                        all_previous_block_begin = all_previous_events(ismember(all_previous_events.value,[30 68]),:);
                        this_block_begin = all_previous_block_begin(end,:);
                        % remove from the output event structure
                        event([event.sample] > this_block_begin.sample & [event.sample] < next_block_begin.sample) = [];
                    end
                case 240
                    % Utilisée seulement pour définir les blocs.
                    % Repair: soit ajouter un trigger artificiel entre les deux triggers voisins (en amont 250 ou 80: fin de l’échelle ou début de session, et en aval 112/150: premier flip du bloc); soit juste remplacer par le premier flip du bloc (ou premier flip du bloc .sample-1 pour ne pas que ce soit exactement les mêmes..?)...
                    all_future_events =  allevents(i_val_allevents_last:end,:);
                    all_future_112_150 = all_future_events(ismember(all_future_events.value,[112 150]),:);
                    first_flip_next_block = all_future_112_150(1,:);
                    event(end+1) = table2struct(first_flip_next_block);
                    event = sortstruct(event,'sample');
            end
        end
    end
    disp('Warning: some triggers missing?!? Attempting a repair...')
    keyboard
    
end

fs = alldata.fsample;

x = [alldata.trial{:}];
