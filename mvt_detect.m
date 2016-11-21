% we get here with
%   mvt_list_trig_act_trials
%   mvt_list_trig_pass_trials
%   data_emg

% we will integrate a filtered & standardized version of data_emg
% between consecutive triggers separated by act & pass, create 2
% distributions, find optimal criterion
cols = colors({'red3','forestgreen'});

% filter and standardize emg
data_emg_filtered = eegfilt(data_emg,fs,20,250);

data_force_detector_filtered = eegfilt(data_force_detector,'srate',fs, 'hicutoff',20); %% low pass filter which blocks frequencies below .5Hz


data_emg2 = data_emg_filtered.^2;
data_emg_z = data_emg2 / std(data_emg2,[],2);


% all triggers
[alltrigs, idx] = sort([mvt_list_trig_act_trials mvt_list_trig_pass_trials]);
alltrigs = round(alltrigs - mvt_det_low_bound);
alltrigs_gps = [ones(1,numel(mvt_list_trig_act_trials)) 2*ones(1,numel(mvt_list_trig_pass_trials))];
alltrigs_gps = alltrigs_gps(idx);

alltrigs(end+1) = alltrigs(end)+round(mean(diff(alltrigs))); % add one 'trigger' at the end
emg_integ = [];force_integ = [];emg_peak = [];
for itri = 1:numel(alltrigs)-1
    idx_loc = [alltrigs(itri):alltrigs(itri+1)] + round(.15*fs);% shift 150ms we consider min RT
    emg_integ(itri) = sum(data_emg_z(idx_loc));
    
    [dum,emg_peak(itri)] = max(data_emg_z(idx_loc));
    emg_peak(itri) = emg_peak(itri) + alltrigs(itri)-1;
    
    force_integ(itri) = sum(data_force_detector_filtered(idx_loc));
end
emg_integ_act = emg_integ(alltrigs_gps == 1);
emg_integ_pass = emg_integ(alltrigs_gps == 2);
force_integ_act = force_integ(alltrigs_gps == 1);
force_integ_pass = force_integ(alltrigs_gps == 2);

[x_mvt,y_mvt,t,auc_mvt,opt_mvt] = perfcurve(alltrigs_gps,emg_integ,1);
best_mvt = t((x_mvt==opt_mvt(1))&(y_mvt==opt_mvt(2)));

mvt_unexp = sum(emg_integ >= best_mvt & alltrigs_gps == 2);% movement when passive trial
mvt_miss = sum(emg_integ < best_mvt & alltrigs_gps == 1);% no move when active trial

% compute RT for correct active trials
mvt_act_corr = emg_integ >= best_mvt & alltrigs_gps == 1;

% mvt RT only for correct active trials
tmp = bsxfun(@minus,emg_peak', alltrigs(mvt_act_corr));
mvt_RT = NaN(1,size(tmp,2));
for i = 1:size(tmp,2)
    mvt_RT(i) = tmp(find(tmp(:,i) > 0,1),i)/fs;
end

[x_force,y_force,t,auc_force,opt_force] = perfcurve(alltrigs_gps,force_integ,1);
best_force = t((x_force==opt_force(1))&(y_force==opt_force(2)));

if doplots
    figure(222)
    set(gcf,'position', [0, 0, 1500, 900])
    subplot(3,1,1);cla
    plot(data_emg_z)
    hold on
    yl = ylim;
    ylim('manual')
    for i = 1:numel(emg_integ)
        text(alltrigs(i),yl(2),num2str(round(emg_integ(i))),'rotation',60)
    end
    % detected conditions based on mvt
    % will only highlight errors
    alltrigs_detected_gps = 2 - [emg_integ>=best_mvt];
    errors = find(alltrigs_gps ~= alltrigs_detected_gps);
    for i = 1:numel(errors)
        fill(alltrigs([errors(i) errors(i) errors(i)+1 errors(i)+1]), repflipped(ylim),'r','FaceAlpha',.2)
    end
    % real conditions
    vline(alltrigs(1:end-1),'color',cols(alltrigs_gps,:),'linewidth',2)
    vline(alltrigs(1:end-1) + round(.15*fs),'color',cols(alltrigs_gps,:),'linewidth',1,'linestyle',':')
    
    xl = xlim;
    text(xl(2),yl(2),['Cutoff ' num2str(round(best_mvt))])
    title({['Active/Passive: ' num2str(expected_nb_active_mvt) '/' num2str(expected_nb_pass_mvt)] 'EMG'})
    
    subplot(3,1,2);cla
    plot(data_force_detector_filtered)
    hold on
    yl = ylim;
    ylim('manual')
    for i = 1:numel(emg_integ)
        text(alltrigs(i),yl(2),num2str(round(force_integ(i)*1e-6)),'rotation',60)
    end
    % detected conditions based on force
    % will only highlight errors
    alltrigs_detected_gps = 2 - [force_integ>=best_force];
    errors = find(alltrigs_gps ~= alltrigs_detected_gps);
    for i = 1:numel(errors)
        fill(alltrigs([errors(i) errors(i) errors(i)+1 errors(i)+1]), repflipped(ylim),'r','FaceAlpha',.2)
    end
    % real conditions
    vline(alltrigs(1:end-1),'color',cols(alltrigs_gps,:),'linewidth',2)
    vline(alltrigs(1:end-1) + round(.15*fs),'color',cols(alltrigs_gps,:),'linewidth',1,'linestyle',':')
    text(xl(2),yl(2),['Cutoff ' num2str(round(best_force*1e-6))])
    title('Force detector')
    
    % show histogram
    bounds = linspace(0,max(emg_integ),20);
    subplot(3,4,[9 10]);cla
    count_emg_integ_act = histc(emg_integ_act,bounds);
    count_emg_integ_pass = histc(emg_integ_pass,bounds);
    [hbars] = bar(bounds+(bounds(2)-bounds(1))/2,[count_emg_integ_act(:) count_emg_integ_pass(:)]);
    xlim([0 max(bounds+(bounds(2)-bounds(1)))])
    hold on
    vline(bounds,':k')
    hthresh = vline(best_mvt,'b');
    yl = ylim;
    scatter(emg_integ_act,repmat(yl(2)-diff(yl)/10,size(emg_integ_act)))
    scatter(emg_integ_act(emg_integ_act<best_mvt),repmat(yl(2)-diff(yl)/10,size(emg_integ_act(emg_integ_act<best_mvt))),[],'r','markerfacecolor','r')
    scatter(emg_integ_pass,repmat(yl(2)-2*diff(yl)/10,size(emg_integ_pass)))
    scatter(emg_integ_pass(emg_integ_pass>=best_mvt),repmat(yl(2)-2*diff(yl)/10,size(emg_integ_pass(emg_integ_pass>=best_mvt))),[],'r','markerfacecolor','r')
    
    legend([hbars hthresh] , {'Active','Passive','Threshold'},'location','best');
    title('Histogram of EMG')
    
    subplot(3,4,11);cla
    scatter(x_mvt,y_mvt);
    hold on
    plot(opt_mvt(1),opt_mvt(2),'or')
    xlabel('False detections')
    ylabel('Correct detections')
    title('ROC curve and optimal criterion')
    
    subplot(3,4,12);cla
    hist(mvt_RT);
    xlabel('EMG peak (RT)')
    title('Histogram of RTs based on EMG')
    waitfor(gcf)
end