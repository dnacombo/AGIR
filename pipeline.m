%% data loading : ONLY IF NEW SET OF DATA
clear
clc;

addpath(fullfile(cd,'general'))
addpath(fullfile(cd,'heart'))
addpath(fullfile(cd,'general','firfilt1.6.1'))
addpath(fullfile(cd,'general', 'plot'))
addpath(fullfile(cd,'fieldtrip'))
ft_defaults

doplots = 0; % whether or not we're going to give visual feedback.
save_each = 0; % whether or not to save one file per block.

% choose yours
% datadir = fullfile('C:\Users\Clémence\Desktop\pipeline_clemence');
datadir = fullfile('/media/max/DATADISK/AGIR_data');

participant = 'participant22';

[fs,x,event,C] = loadallgdf(datadir,participant);




%% STEP 1//// identifying blocks and lists of triggers
%
trigger_list = [event.value]; %% to have real value of the list of sent triggers
sample_trigger = [event.sample];

% columns in the matrix "C".
col_nb_total_trial=1;
col_sync=2;
col_percentage = 3;
col_this_trial_type = 4;
col_what_proportion_tells_type = 5;

%specific triggers lists in the data set
% detected HB
list_trig_det_HB=[event(trigger_list == 200 ).sample];
% visual warning (Actif/Passif
list_trig_flips=[event(trigger_list == 112 | trigger_list == 150 ).sample];
% Trials actifs
list_trig_act_trials =[event(trigger_list == 150 ).sample];
% Trials passifs
list_trig_pass_trials =[event(trigger_list == 112 ).sample];
% 1st trial begins block
list_trig_first_trials = [event(trigger_list == 240).sample];
% defining synchronicity for each block (end of block except for last one).
list_trig_synty_def=[event(trigger_list == 30 | trigger_list == 68 ).sample];
% all triggers for block definition
list_trig_for_blocks = sort([list_trig_synty_def, list_trig_first_trials]);
% passive mvt done
list_trig_passive_mvt_done = [event(trigger_list == 130).sample];

%refine the triggers list for identifying the blocks within each dataset

% Max: I comment the ifs below because we have all triggers now (all datasets at
% once) so we need to remove 1st trigger of block 1 and add the last event
% as if this_sata_set == 7

% if this_data_set == 1
list_trig_for_blocks(1) = []; % si c'est le premier gdf vire le 1er trigger (qui est bloc def)
% elseif this_data_set == 7
list_trig_for_blocks = [list_trig_for_blocks, [event(trigger_list == 220 ).sample]]; %% careful, last block will be too long..
% clemence = rajoute le trigger de fin de la t�che si c'est le
% dernier bloc
% end
%checker dans le code de l'exp� mais il semble que dans tous les autres gdf il manque une borne

% cle_first_trials = [trigger_list(trigger_list == 240) ; event(trigger_list == 240).sample]
% cle_blocks_def = [trigger_list(trigger_list == 30|trigger_list==68) ; event(trigger_list == 30|trigger_list==68).sample]
%set = find(trigger_list==240 |trigger_list==230 |trigger_list==30 |trigger_list==68 | trigger_list==80);
%list = [set; trigger_list(set)]'
[blocks_def] = identifying_blocks(list_trig_for_blocks);

BigTable = struct('suj',{});
for this_block = 1:size(blocks_def,1)
    
    
    %% STEP 2 ///// this block properties
    
    filename = [participant '_block_',num2str(this_block)];
    
    %what mvt are expected
    if C{this_block}(1, col_what_proportion_tells_type) == 1
        expected_nb_active_mvt = C{this_block}(1, col_nb_total_trial)*10*C{this_block}(1, col_percentage);
        expected_nb_pass_mvt = C{this_block}(1, col_nb_total_trial)*10 - expected_nb_active_mvt;
    elseif C{this_block}(1, col_what_proportion_tells_type) == 2
        expected_nb_active_mvt = C{this_block}(1, col_nb_total_trial)*10*(1-C{this_block}(1, col_percentage));
        expected_nb_pass_mvt = C{this_block}(1, col_nb_total_trial)*10 - expected_nb_active_mvt;
    end
    
    % to plot a general figure of the block
    
    if doplots
        figure(1)
        
        timetoplot = blocks_def(this_block, 1 ) - fs/2 : blocks_def(this_block, 2);
        toplot = x(:,timetoplot);
        toplot = baseline_corr(toplot,2,1);
        
        subplot(4,1,1);
        plot(timetoplot,toplot(1,:))
        title('Trigger channel')
        subplot(4,1,2);
        plot(timetoplot,toplot(2,:))
        title('ECG')
        subplot(4,1,3);
        plot(timetoplot,toplot(4,:)-toplot(3,:))
        title('EMG')
        subplot(4,1,4);
        plot(timetoplot,toplot(5,:))
        title('Force')
        drawnow
        
        xlabel('Close this window if everything is ok.')
        waitfor(gcf)
    end
    %save(filename);
    
    
    
    %% STEP 3 //// Raw HB signal analysis, HB detection from raw signal
    
    % block definition : take a large data which will be usefull for all
    % later analysis with HB. We want to be sure to get the HB before and
    % after the block as well
    % 7 + 2s before the blocks begin, for asycnrhoneous conditions
    % if this brings us before first sample we use 1
    low_bound = max(1,blocks_def(this_block, 1 ) - 2*fs - 7*fs); 
    high_bound =  blocks_def(this_block, 2) + 2*fs; %% 1s after the block ends, for later normalization of intervals according to IBI
    data_sample_heart = x(2,[low_bound : high_bound]);
    
    % HB peak detection
    cfg = [];
    cfg.downsample = 'no';
    cfg.downrate = 300;
    cfg.fsample = fs;
    cfg.hplocut = 1;
    cfg.plotthresh = doplots;
    cfg.plotcorr = doplots;
    
    [HB_det_raw_signal] = heart_peak_detect(cfg,data_sample_heart); %%  R peaks in this data sample
    HB_det_raw_signal = [HB_det_raw_signal.R_sample];
    % [HB_det_raw_signal_old] = heart_peak_detection_old2(data_sample_heart,fs); %%  R peaks in this data sample
    % below is done within heart_peak_detect
    % HB_det_raw_signal_old = round(HB_det_raw_signal_old*fs/300); % resample with original sampling
    
    HB_det_raw_signal = HB_det_raw_signal +  low_bound; % real values of samples
    
    % general IBI in the block
    nuIBI_intervals = diff(HB_det_raw_signal);
    IBI_mean = mean(nuIBI_intervals);
    
    %             close all;
    
    %save(filename);
    
    
    
    
    %% STEP 4 ////// efficiency of HB on-line detection by fieldtrip
    
    % block definition for efficiency of HB detection
    if C{this_block}(1, col_sync) ==2 %% if Async, look at HB detection only for time window of the block - 7s
        eff_low_bound = blocks_def(this_block, 1 ) -7*fs - .5*fs ; %% the .5*fs is to get the first HB
        eff_high_bound = blocks_def(this_block, 2 ) -7*fs ; %% 7 seconds for the delay in async and 1 sec for the big window I took to get the last HB for normalizing afterwards
    elseif C{this_block}(1, col_sync) ==1 %% if Sync, look at HB detection for time window of the block
        eff_low_bound = blocks_def(this_block, 1) - .5*fs ;
        eff_high_bound = blocks_def(this_block, 2 );
    end
    
    % we assume all HB were detected offline.
    % How many were also detected online?
    % calculating the efficiency = the % of HB that were detected
    eff_list_trig_det_HB = list_trig_det_HB(list_trig_det_HB >= eff_low_bound & list_trig_det_HB <= eff_high_bound); % triggers detecting HB
    eff_list_HB_det_raw_signal = HB_det_raw_signal(HB_det_raw_signal>=eff_low_bound & HB_det_raw_signal<=eff_high_bound); % real HB
    
    [eff_online2offline, eff_HB_detection_times] = dsearchn(eff_list_HB_det_raw_signal', eff_list_trig_det_HB');
    eff_proportion_det_HB = 100*numel(eff_online2offline) / numel(eff_list_HB_det_raw_signal);
    eff_nb_of_missed_HB = numel(eff_list_HB_det_raw_signal) - numel(eff_list_trig_det_HB);
    eff_mean_det_time_sample = mean(eff_HB_detection_times);
    eff_mean_det_time_secondes = eff_mean_det_time_sample/fs;
    
    %             eff_proportion_det_HB = 100 - (length( eff_list_HB_det_raw_signal) - length( eff_list_trig_det_HB))*100/length( eff_list_HB_det_raw_signal);
    %             eff_nb_of_missed_HB =  length( eff_list_HB_det_raw_signal) - length( eff_list_trig_det_HB);
    %
    %             %calculating the mean of the detection time for the detected HB
    %             eff_HB_detection_times = [];
    %             for j = 1 : length(eff_list_trig_det_HB)
    %                 r = abs(eff_list_trig_det_HB(j) - eff_list_HB_det_raw_signal);
    %
    %                 if min(r) < .35*fs %% assuming all detection occur before 350ms
    %                     eff_HB_detection_times = [eff_HB_detection_times, min(r)];
    %                 end
    %             end
    %             eff_mean_det_time_sample = mean(eff_HB_detection_times);
    %             eff_mean_det_time_secondes = eff_mean_det_time_sample/fs;
    %save(filename);
    
    
    
    %% STEP 5 ////// testing synchrony vs asynchrony of flips and heartbeats.
    
    
    % try
    % block definition for testing synchrony vs asynchrony of heartbeats.
    Synty_low_bound = blocks_def(this_block, 1 ) - 2*fs;
    Synty_high_bound = blocks_def(this_block, 2 ) + 2*fs; %% this is to be sure to get 1 HB after the last mvt, for normalization afterwards.
    %but this should not apply to list_trig_flips or we start getting flips for
    %the next/previous block...
    
    Synty_list_trig_flips = list_trig_flips(list_trig_flips >= Synty_low_bound & list_trig_flips <= Synty_high_bound);
    Synty_HB_det_raw_signal = HB_det_raw_signal(HB_det_raw_signal >= Synty_low_bound & HB_det_raw_signal <= Synty_high_bound);
    
    Synty_delays_flip_closest_last_HB = [];
    Synty_delays_flip_closest_next_HB = [];
    normalized_delay_flip_to_HB = [];
    
    
    % this is to guess the position of flip triggers if missing, using the
    % timing of flips in the C matrix. For a wiring reason, all flips
    % triggers were not sent though flips were desplayed on the screen.
    
    if length(Synty_list_trig_flips) < 10*(C{this_block}(1, col_nb_total_trial)) %% if, somehow, all flips triggers have not been sent
        for k = 1 : 10*(C{this_block}(1, col_nb_total_trial)) -1
            a = C{this_block}(k + 1, 6)*fs - C{this_block}(k, 6)*fs;
            if k+1 <= length(Synty_list_trig_flips)
                b = Synty_list_trig_flips(k+1) - Synty_list_trig_flips(k);
                if b > a + .25*fs   %% if the next flip is within a window of 500ms around the expected time
                    c = Synty_list_trig_flips(k) + a;
                    Synty_list_trig_flips = [Synty_list_trig_flips, c];
                end
                Synty_list_trig_flips = sort(Synty_list_trig_flips);
            else
                c = Synty_list_trig_flips(k) + a;
                Synty_list_trig_flips = [Synty_list_trig_flips, c];
            end
        end
    end
    
    
    % now that we are sure to have the good list of flip triggers,
    % normalize the delay within the IBI of each trial.
    
    % normalized HB-flip-HB
    tmp = bsxfun(@minus,Synty_list_trig_flips', HB_det_raw_signal);
    previousHB = NaN(size(tmp,1),1);
    nextHB = NaN(size(tmp,1),1);
    for i = 1:size(tmp,1)
        try
            previousHB(i) = tmp(i,find(tmp(i,:) > 0,1,'last'))/fs;
            nextHB(i) = tmp(i,find(tmp(i,:) < 0,1,'first'))/fs;
        end
    end
    normalized_delay_flip_to_HB = previousHB' ./ (previousHB-nextHB)';
    %
    % for j = 1 : length(Synty_list_trig_flips)
    %     r = Synty_list_trig_flips(j) - Synty_HB_det_raw_signal;
    %     r = r(r > 0);
    %     Synty_delays_flip_closest_last_HB = [Synty_delays_flip_closest_last_HB, min(r)];
    %     z = Synty_list_trig_flips(j) - Synty_HB_det_raw_signal;
    %     z = z(z < 0);
    %     z =abs(z);
    %     Synty_delays_flip_closest_next_HB = [Synty_delays_flip_closest_next_HB, min(z)];
    %
    %     if not(isempty(z))% if there's no
    %         normalized_delay_flip_to_HB = [normalized_delay_flip_to_HB, min(r)/(min(r)+min(z))];
    %     end
    %
    % end
    
    
    %now compare the normalized data to a continuous
    %distribution : create a histogram of 10 bins of the normalized data
    %and then compare it to the same histogram for continuous
    %distribution
    
    %     bin_normalized_data(1, 1:10) =[0];
    %     for k = 1 : length(normalized_delay_flip_to_HB)
    %         for m =  1 : 10
    %             if  normalized_delay_flip_to_HB(k) >(m-1)/10 & normalized_delay_flip_to_HB(k) < m/10
    %                 bin_normalized_data(1,m) = [bin_normalized_data(1,m) + 1];
    %             end
    %         end
    %     end
    binranges = 0:.1:1;
    bin_normalized_data = histc(normalized_delay_flip_to_HB,binranges);
    bin_normalized_data(end-1) = bin_normalized_data(end-1)+bin_normalized_data(end);
    bin_normalized_data(end) = [];
    
    % expected if distribution perfectly flat (given how many flip_to_HB
    % measures we have rather than how many we should have)
    continuous_distribution_trials_per_bin = sum(bin_normalized_data) / numel(bin_normalized_data);
    bin_normalized_continuous = repmat(continuous_distribution_trials_per_bin,1,numel(bin_normalized_data));
    
    distribution_difference = bin_normalized_data - bin_normalized_continuous;
    
    continuous_distribution_out_trials = sum(distribution_difference(distribution_difference>0));
    
    % continuous_distribution_out_trials = 0;
    % continuous_distribution_trials_per_bin = C{this_block_in_general}(1, col_nb_total_trial);
    % for k = 1 : length(bin_normalized_data)
    %     if bin_normalized_data(1, k) > continuous_distribution_trials_per_bin %% if the value in this bin of the normalized data > the expected value if the repartition were continuous
    %         continuous_distribution_out_trials = continuous_distribution_out_trials + (bin_normalized_data(1, k) - continuous_distribution_trials_per_bin);
    %     end
    % end
    
    %score for global syncrhonicity, normalized to 0-100%
    synchronicity_HB_flips = (continuous_distribution_out_trials*100 / (C{this_block}(1, col_nb_total_trial)*10 - continuous_distribution_trials_per_bin)); %% gives the proportion of good disparity, max 100% but min > 0%
    
    
    %save(filename);
    
    
    % figure
    % hist(Synty_delays_flip_closest_last_HB);
    
    %             figure
    %             plot(normalized_delay_flip_to_HB, 'ob');
    %             axis([0 length(normalized_delay_flip_to_HB) 0 1])
    %
    %             title(['synchronicity score : ', num2str(synchronicity_HB_flips)]);
    
    
    %% STEP 6 ///// Active mvt detection
    
    % try
    % block data sample definition
    mvt_det_low_bound = Synty_low_bound;
    mvt_det_high_bound =  Synty_high_bound;
    data_force_detector = x(5,[mvt_det_low_bound : mvt_det_high_bound]);
    data_emg = x(4,[mvt_det_low_bound : mvt_det_high_bound])- x(3,[mvt_det_low_bound : mvt_det_high_bound]);
    
    mvt_list_trig_act_trials = [];
    mvt_list_trig_pass_trials = [];
    mvt_list_trig_flips= Synty_list_trig_flips;   %% visual warning flips samples, corrected in needed by the flips timing
    
    %separate flips active and passive in the flip triggers list
    for k = 1 : C{this_block}(1, col_nb_total_trial)*10
        if C{this_block}(1, col_what_proportion_tells_type) == 1
            if C{this_block}(k, col_this_trial_type) ==1
                mvt_list_trig_act_trials = [mvt_list_trig_act_trials, mvt_list_trig_flips(k)];
            elseif C{this_block}(k, col_this_trial_type) ==2
                mvt_list_trig_pass_trials = [mvt_list_trig_pass_trials, mvt_list_trig_flips(k)];
            end
        elseif C{this_block}(1, col_what_proportion_tells_type) == 2
            if C{this_block}(k, col_this_trial_type) ==2
                mvt_list_trig_act_trials = [mvt_list_trig_act_trials, mvt_list_trig_flips(k)];
            elseif C{this_block}(k, col_this_trial_type) ==1
                mvt_list_trig_pass_trials = [mvt_list_trig_pass_trials, mvt_list_trig_flips(k)];
            end
        end
    end
    
    %function for mvt cetection with EMG only.
    mvt_detect
    
    
    %save(filename);
    
    
    
    
    %% STEP 7 ///// synchronicity of mvt and motor with heartbeats.
    
    % normalized HB-EMGpeak-HB
    allmvt = emg_peak;
    % remove missed movements
    allmvt([alltrigs_gps == 1 & emg_integ < best_mvt]) = NaN;
    
    %             allmvt(alltrigs_gps==2) = list_trig_passive_mvt_done(list_trig_passive_mvt_done >= Synty_low_bound & list_trig_passive_mvt_done <= Synty_high_bound)  - mvt_det_low_bound;
    pass_run = find(alltrigs_gps == 2);
    flip = mvt_list_trig_pass_trials  - mvt_det_low_bound;
    trigs_mvt_done_loc = list_trig_passive_mvt_done' - mvt_det_low_bound;
    for i = 1:numel(flip)
        [closest,distance] = dsearchn(trigs_mvt_done_loc,flip(i));
        if distance < .7 * fs
            allmvt(pass_run(i)) = trigs_mvt_done_loc(closest);
        else
            allmvt(pass_run(i)) = NaN;
        end
    end
    
    tmp = bsxfun(@minus,allmvt', HB_det_raw_signal - mvt_det_low_bound);
    previousHB = NaN(size(tmp,1),1);
    nextHB = NaN(size(tmp,1),1);
    for i = 1:size(tmp,1)
        try
            previousHB(i) = tmp(i,find(tmp(i,:) > 0,1,'last'))/fs;
            nextHB(i) = tmp(i,find(tmp(i,:) < 0,1,'first'))/fs;
        end
    end
    HB2allmvt = previousHB ./ (previousHB-nextHB);
    cols = colors({'red','grey30'});
    
    
    
    
    binranges = 0:.1:1;
    bin_normalized_data = histc(HB2allmvt',binranges);
    bin_normalized_data(end-1) = bin_normalized_data(end-1)+bin_normalized_data(end);
    bin_normalized_data(end) = [];
    
    continuous_distribution_out_trials = 0;
    continuous_distribution_trials_per_bin = C{this_block}(1, col_nb_total_trial);
    for k = 1 : length(bin_normalized_data)
        if bin_normalized_data(1, k) > continuous_distribution_trials_per_bin %% if the value in this bin of the normalized data > the expected value if the repartition were continuous
            continuous_distribution_out_trials = continuous_distribution_out_trials + (bin_normalized_data(1, k) - continuous_distribution_trials_per_bin);
        end
    end
    
    
    %score for global syncrhonicity, normalized to 0-100%
    synchronicity_HB_mvt = (continuous_distribution_out_trials*100 / (C{this_block}(1, col_nb_total_trial)*10 - continuous_distribution_trials_per_bin)); %% gives the proportion of good disparity, max 100% but min > 0%
    
    
    if doplots
        figure(444)
        subplot(121)
        plot([normalized_delay_flip_to_HB],'o');
        ylim([0 1]);
        title('Normalized HB 2 flip')
        subplot(122)
        h = scatter(1:numel(HB2allmvt),[HB2allmvt],40,cols(alltrigs_gps,:)); 
        ylim([0 1]);
        title({'Normalized HB 2 movement' 'red is active'})
        waitfor(gcf)
    end
    
    
    % cle: adding Victor's flow_flips
    list_IFI = [];
    for i = 1 : size(C{this_block},1) - 1
        list_IFI = [list_IFI, C{this_block}(i+1, 6) - C{this_block}(i, 6) ];
    end
    nb_HBflip_skipped = sum(list_IFI>1.33);%% considering if the interval between 2 flips exceeds 1.3 s, then at least one HB was skipped
    % cle: adding sum() so we get the nb.
    
    BigTable(end+1).suj = participant;
    BigTable(end).block = this_block;
    BigTable(end).HB_missed_online = eff_nb_of_missed_HB;
    BigTable(end).HB_error_prop = eff_proportion_det_HB;
    BigTable(end).nb_HBflip_skipped = nb_HBflip_skipped;
    BigTable(end).mvt_missed = mvt_miss;
    BigTable(end).mvt_unexp = mvt_unexp;
    BigTable(end).synch_HB_flips = synchronicity_HB_flips;
    BigTable(end).synch_HB_mvt = synchronicity_HB_mvt;
    BigTable(end).sync_def_was_ok = C{this_block}(1,9);
    
    if save_each
        vars = who;
        vars(regexpcell(vars,'data_.*|x|BigTable')) = [];
        save(filename,vars{:});
    end
    
end

save(participant,'BigTable')

