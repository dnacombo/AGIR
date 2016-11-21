function [DATA STUDY ALLEEG] = std_get_readdata_tmp(STUDY,ALLEEG,varargin)

% [DATA STUDY ALLEEG] = std_get_readdata(STUDY,ALLEEG,varargin)
%
% go retrieve single trial data already read from STUDY (i.e. stored inside
% the structure) uses the same syntax as std_readersp
% note that 'singletrials' must be 'on'
%
% returned data is a {case x var1 x var2} cell array containing one cell per
% dataset, itself containing the AVERAGED data on a given time-frequency-channel_or_component window.
% trials x cluster_or_changrp
%
% if 'inject' is true, then data is injected in ALLEEG.event
% one data point per trial in injectname field
%

STUDY = pop_erspparams( STUDY, 'default');
STUDY = pop_erpimparams(STUDY, 'default');
[opt] = finputcheck( varargin, { ...
    'inject'        'boolean' []             1
    'injectname'    'string'  []            'injected_data' ;
    'design'        'integer' []             STUDY.currentdesign;
    'channels'      'cell'    {}             {};
    'clusters'      'integer' []             [];
    'freqrange'     'real'    []             [0 inf];
    'timerange'     'real'    []             [-inf inf];
    'singletrials'  'string'  { 'on' }       'on';
    'component'     'integer' []             [];
    'datatype'      'string'  { 'ersp'}     'ersp'; ...,'itc','erpim' 'timef'
    'subject'       'string'  []             '' }    );

if isstr(opt), error(opt); end;
if strcmpi(opt.datatype, 'erpim'), if isnan(opt.timerange), opt.timerange = STUDY.etc.erpimparams.timerange; end;
    opt.freqrange = opt.trialrange;
    ordinate      = 'trials';
else                               if isnan(opt.timerange) opt.timerange = STUDY.etc.erspparams.timerange; end;
    ordinate      = 'freqs';
end;
nc = max(length(STUDY.design(opt.design).variable(1).value),1);
ng = max(length(STUDY.design(opt.design).variable(2).value),1);
paired1 = STUDY.design(opt.design).variable(1).pairing;
paired2 = STUDY.design(opt.design).variable(2).pairing;

dtype = opt.datatype;

if nargout < 2
    opt.inject = 0;
    % no use if we're not asking for results
end

if isempty(opt.clusters)
    % if we're not asking for clusters, we'll lookfor channels
    clustidx = chnb(opt.channels);
    changrp_or_clust = 'changrp';
else
    % looking for clusters
    clustidx = opt.clusters;
    changrp_or_clust = 'cluster';
end

desidx = opt.design;
subjidx = opt.subject;
cases = STUDY.design(desidx).cases.value;
if isempty(subjidx)
    subjidx = cases;
end
var1 = STUDY.design(desidx).variable(1).value;
var2 = STUDY.design(desidx).variable(2).value;

% trials all datasets
sets_nb_trials = [ALLEEG.trials];
if opt.inject
    for i = 1:numel(ALLEEG)
        fs = fieldnames(ALLEEG(i).event);
        ifs = strncmp(opt.injectname,fs,numel(opt.injectname));
        ALLEEG(i).event = rmfield(ALLEEG(i).event,fs(ifs));
    end
end

DATA = cell(numel(cases),numel(var1),numel(var2));
for i = 1:numel(DATA)
    DATA{i} = cell(1,numel(clustidx));
end

% trialdat will reorganize the data so that we can in the end
% target specific datasets & variables
% for each value of var1
for ivar1= 1:numel(var1)
    % for each value of var2
    for ivar2 = 1:numel(var2)
        for i_clust = 1:numel(clustidx)
            timeidx = timepts(opt.timerange, STUDY.(changrp_or_clust)(clustidx(i_clust)).([dtype 'times']));
            freqidx = timepts(opt.freqrange,STUDY.(changrp_or_clust)(clustidx(i_clust)).([dtype 'freqs']));
            % we pick each component
            icomploc = zeros(1,numel(STUDY.datasetinfo));
            for icomp = 1:size(STUDY.(changrp_or_clust)(clustidx(i_clust)).([dtype 'subjinds']){ivar1,ivar2},2)
                cellidx = STUDY.(changrp_or_clust)(clustidx(i_clust)).setinds{ivar1,ivar2}(icomp);
                setidx = STUDY.design(desidx).cell(cellidx).dataset;
                caseidx = strcmp(STUDY.design(desidx).cases.value, ...
                    STUDY.design(desidx).cell(cellidx).case) & ...
                    strcmp(STUDY.design(desidx).cases.value, ...
                    subjidx);
                if all(~caseidx)
                    continue
                end
                tricompidx = STUDY.(changrp_or_clust)(clustidx(i_clust)).([dtype 'subjinds']){ivar1,ivar2}(:,icomp);
                tricompidx = tricompidx(1):tricompidx(2);
                setsthiscell = unique([STUDY.design(desidx).cell(cellidx).dataset]);
                clear trialorder
                for iset = 1:numel(setidx)
                    icomploc(setidx(iset)) = icomploc(setidx(iset)) + 1;
                    if isempty(DATA{caseidx,ivar1,ivar2}{i_clust})
                        for i = 1:numel(setidx)
                            DATA{caseidx,ivar1,ivar2}{i_clust}{i} = NaN(sum( [STUDY.design(desidx).cell([STUDY.changrp(clustidx(i_clust)).setinds{ivar1,ivar2}]).dataset] == setidx(iset)),sets_nb_trials(setidx(iset)),'single');
                        end
                    end
%                     thissetthiscell = find(setidx(iset) == setsthiscell);
%                     icompcell = icomp - sum(ismember([STUDY.design(desidx).cell(cellidx).dataset(iset)],setsthiscell(1:thissetthiscell-1)));
                    trialorderthiscell{iset} = STUDY.design(desidx).cell(cellidx).trials{iset};
                    itrithiscellthisset = (1:numel(trialorderthiscell{iset})) + numel([trialorderthiscell{1:iset-1}]);
                    % we gather the data here for specific cases and
                    % sets within var1,var2
                    DATA{caseidx,ivar1,ivar2}{i_clust}{iset}(icomploc(setidx(iset)),trialorderthiscell{iset}) = permute(mean(mean(STUDY.(changrp_or_clust)(clustidx(i_clust)).([dtype 'datatrials']){ivar1,ivar2}(timeidx,freqidx,tricompidx(itrithiscellthisset)),1),2),[3 1 2]);
                    if opt.inject && icomploc(setidx(iset)) == sum([STUDY.design(desidx).cell([STUDY.changrp(clustidx(i_clust)).setinds{ivar1,ivar2}]).dataset] == setidx(iset))
                        % average the component/chan dimension
                        tmpdat = nanmean(DATA{caseidx,ivar1,ivar2}{i_clust}{iset},1);
                        
                        
                        epochs = [ALLEEG(setidx(iset)).event.epoch];
                        tmpdatallep = tmpdat(epochs);
                        switch changrp_or_clust
                            case 'changrp'
                                [~,~, ch] = chnb(clustidx(i_clust));
                            case 'cluster'
                                ch = ['cluster' num2str(clustidx(i_clust),'%02d')];
                        end
                        injfield = [opt.injectname '_' ch];
                        if nargout == 3
                            for i_ev = find(~isnan(tmpdatallep(:)))'
                                ALLEEG(setidx(iset)).event(i_ev).(injfield) = tmpdatallep(i_ev);
                            end
                        end
                        for i_tri = find(~isnan(tmpdat))
                            STUDY.datasetinfo(setidx(iset)).trialinfo(i_tri).(injfield) = tmpdat(i_tri);
                        end
                    end
%                     trialsetnum(caseidx,ivar1,ivar2) = setidx(iset);
                end
            end
        end
    end
end


 