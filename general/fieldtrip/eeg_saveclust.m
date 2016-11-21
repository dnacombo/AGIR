function eeg_saveclust(STUDY,clustname,del)

% eeg_saveclust(STUDY,clustname,del)
% export current clusters from STUDY so you can play around a bit.
% overwrite any previously saved cluster with same name.
% delete cluster of that name if del is true.
% delete clusters file if del is true and clustname is empty.
%

if not(exist('del','var')) || ~del
    del = 0;
else
    del = 1;
end
studyfilename = fullfile(STUDY.filepath,STUDY.filename);
filename = strrep(studyfilename,'.study','.clusters');
% we don't want to save the data so we clear all big fields.
% might have to add some more here.
toremove = {
    'erspfreqs'
    'erspdatatrials'
    'erspsubjinds'
    'ersptimes'
    'ersptrialinfo'};
for i = 1:numel(toremove)
    try
        STUDY.cluster = rmfield(STUDY.cluster,toremove);
    end
end

try
    load(filename,'-mat');
catch
    CLUSTERS = {};
    clusternames = {};
end
id = strcmp(clusternames,clustname);
if del
    if isempty(clustname)
        delete(filename);
        return
    else
        CLUSTERS(id) = [];
        clusternames(id) = [];
    end
else
    
    if all(~id)
        id = numel(CLUSTERS)+1;
    end
    CLUSTERS{id} = STUDY.cluster;
    clusternames{id} = clustname;
end

save(filename,'CLUSTERS','clusternames');

disp('clusters saved.')

disp('Saving dipfits as well...')
std_savedipfit(STUDY,clustname,del);


