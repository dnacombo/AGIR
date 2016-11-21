function data = std_get_readdata(STUDY,ALLEEG,varargin)



STUDY = pop_erspparams( STUDY, 'default');
STUDY = pop_erpimparams(STUDY, 'default');
[opt moreopts] = finputcheck( varargin, { ...
    'design'        'integer' []             STUDY.currentdesign;
    'channels'      'cell'    {}             {};
    'clusters'      'integer' []             [];
    'freqrange'     'real'    []             [0 inf];
    'timerange'     'real'    []             [-inf inf];
    'singletrials'  'string'  { 'on','off' } 'off';
    'subbaseline'   'string'  { 'on','off' } STUDY.etc.erspparams.subbaseline;
    'component'     'integer' []             [];
    'datatype'      'string'  { 'ersp'} 'ersp'; ...,'itc','erpim' 'timef'
    'subject'       'string'  []             '' }, ...
    'std_readersp', 'ignore');

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

if isempty(opt.clusters)
    % if we're not asking for clusters, we'll lookfor channels
    clustidx = chnb(opt.channels);
    comps_or_clust = 'changrp';
else
    % looking for clusters
    clustidx = opt.clusters;
    comps_or_clust = 'cluster';
end
desidx = opt.design;
subjidx = opt.subject;
cases = STUDY.design(desidx).cases.value;
if isempty(subjidx)
    subjidx = cases;
end
var1 = STUDY.design(desidx).variable(1).value;
var2 = STUDY.design(desidx).variable(2).value;

timeidx = timepts(opt.timerange, STUDY.(comps_or_clust)(clustidx(1)).([dtype 'times']));
freqidx = timepts(opt.freqrange,STUDY.(comps_or_clust)(clustidx(1)).([dtype 'freqs']));

trialdat = cell(numel(cases),numel(var1),numel(var2));
% trialdat will reorganize the data so that we can in the end
% target specific datasets & variables
% for each value of var1
for ivar1= 1:numel(var1)
    % for each value of var2
    for ivar2 = 1:numel(var2)
        % we pick each component
        for icomp = 1:size(STUDY.(comps_or_clust)(clustidx).erspsubjinds{ivar1,ivar2},2)
            cellidx = STUDY.(comps_or_clust)(clustidx).setinds{ivar1,ivar2}(icomp);
            setidx = STUDY.design(desidx).cell(cellidx).dataset;
            caseidx = strcmp(STUDY.design(desidx).cases.value, ...
                STUDY.design(desidx).cell(cellidx).case) & ...
                strcmp(STUDY.design(desidx).cases.value, ...
                subjidx);
            if all(~caseidx)
                continue
            end
            tricompidx = STUDY.(comps_or_clust)(clustidx).erspsubjinds{ivar1,ivar2}(:,icomp);
            tricompidx = tricompidx(1):tricompidx(2);
            setshere = unique([STUDY.design(desidx).cell(cellidx).dataset]);
            clear trialorder
            for iset = 1:numel(setidx)
                thissethere = find(setidx(iset) == setshere);
                icomploc = icomp - sum(ismember([STUDY.design(desidx).cell(cellidx).dataset(iset)],setshere(1:thissethere-1)));
                trialorder{iset} = STUDY.design(desidx).cell(cellidx).trials{iset};
                itriloc = (1:numel(trialorder{iset})) + numel([trialorder{1:iset-1}]);
                % we gather the data here for specific cases and
                % sets within var1 (= a group of conditions)
                trialdat{caseidx,ivar1,ivar2}{iset}(:,:,icomploc,trialorder{iset}) = STUDY.(comps_or_clust)(clustidx).([dtype 'datatrials']){ivar1,ivar2}(timeidx,freqidx,tricompidx(itriloc));
                trialsetnum{caseidx,ivar1,ivar2}(iset) = setidx(iset);
            end
        end
    end
end

data = trialdat;