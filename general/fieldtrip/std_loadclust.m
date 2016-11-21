function varargout = std_loadclust(STUDY, ALLEEG,clustname)

% [STUDY ALLEEG] = std_loadclust(STUDY,ALLEEG, clustname)
% load cluster information stored with tag clustname
% will also load the corresponding dipfit information in all EEG sets.
%
% eeg_loadclust(STUDY) will list possible stored clustname values.


filename = fullfile(STUDY.filepath,strrep(STUDY.filename,'.study','.clusters'));

try
    load(filename,'-mat');
catch
    error('No clusters previously saved')
end
if not(exist('clustname','var'))
    for i = 1:numel(clusternames)
        disp(clusternames{i})
    end
    varargout = {[]};
    return
end
id = strcmp(clusternames,clustname);
if all(~id)
    error('Unknown cluster name')
end
STUDY.cluster = [];
STUDY.cluster = CLUSTERS{id};
disp('Clusters loaded')

disp('Do you also want to (R)eload previously saved dipole fits?')
disp('or do you want to (S)kip that step and keep the dipfits you have?')
s = input('R/S ? ','s');
switch lower(s)
    case 's'
        disp('Skip loading dipole fits. No dataset file has been modified.')
    case 'r'
        disp('Reloading saved dipole fits.')
        [STUDY ALLEEG] = std_loaddipfit(STUDY,ALLEEG,clustname);
end

varargout{1} = STUDY;
varargout{2} = ALLEEG;

