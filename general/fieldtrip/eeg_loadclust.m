function varargout = eeg_loadclust(STUDY,clustname)

filename = fullfile(STUDY.filepath,strrep(STUDY.filename,'.study','.clusters'));

try
    load(filename,'-mat');
catch
    error('No clusters previously saved')
end
if not(exist('clustname','var'))
    disp(clusternames)
    varargout = {[]};
    return
end
id = strcmp(clusternames,clustname);
if all(~id)
    error('Unknown cluster name')
end
STUDY.cluster = [];
STUDY.cluster = CLUSTERS{id};

disp('Loading dipfits as well...')
std_loaddipfit(STUDY,clustname,del);

varargout{1} = STUDY;

