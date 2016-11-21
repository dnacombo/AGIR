function [STUDY ALLEEG] = std_loaddipfit(STUDY,ALLEEG,fitname)

% [STUDY ALLEEG] = std_loaddipfit(STUDY,ALLEEG,fitname)
% load dipfit data stored with tag fitname for all the datasets in this
% study. 
% We actually rewrite all datasets with the loaded dipfit info. this is a bit
% risky...



filename = fullfile(STUDY.filepath,strrep(STUDY.filename,'.study','.dipfits'));

try
    load(filename,'-mat');
catch
    error('No dipfit previously saved')
end
if not(exist('fitname','var'))
    disp(dipfitnames)
    return
end
id = strcmp(dipfitnames,fitname);
if all(~id)
    error('Unknown cluster name')
end
disp('I''ll now write the dipfit information you saved previously to all datasets of the study.')
s = input('You sure you want to proceed? (y/n) ','s');
if ~strcmp(s,'y')
    disp('Abort. No dataset file has been modified.')
    return
end
    
for i_dat = 1:numel(STUDY.datasetinfo)
    EEG = pop_loadset('filepath',STUDY.datasetinfo(i_dat).filepath,'filename',STUDY.datasetinfo(i_dat).filename,'loadmode','info');
    
    EEG.dipfit = DIPFITS{id}{i_dat};
    pop_saveset(EEG,'savemode','resave');
    ALLEEG(i_dat) = EEG;
end


