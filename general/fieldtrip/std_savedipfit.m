function std_savedipfit(STUDY,fitname,del)

% std_savedipfit(STUDY,fitname,del)
% export current dipfits from STUDY so you can play around a bit.
% overwrite any previously saved dipfit with same name.
% delete dipfit of that name if del is true.
% delete whole dipfit file if del is true and fitname is empty.
%


if not(exist('del','var')) || ~del
    del = 0;
else
    del = 1;
end
studyfilename = fullfile(STUDY.filepath,STUDY.filename);
filename = strrep(studyfilename,'.study','.dipfits');
% we don't want to save the data so we clear all big fields.
% might have to add some more here.
try
    load(filename,'-mat');
catch
    DIPFITS = {};
    dipfitnames = {};
end
id = strcmp(dipfitnames,fitname);
if del
    if isempty(fitname)
        delete(filename);
        return
    else
        DIPFITS(id) = [];
        dipfitnames(id) = [];
    end
else
    
    if all(~id)
        id = numel(DIPFITS)+1;
    end
    for i_dat = 1:numel(STUDY.datasetinfo)
        EEG = pop_loadset('filepath',STUDY.datasetinfo(i_dat).filepath,'filename',STUDY.datasetinfo(i_dat).filename,'loadmode','info');
    
        DIPFITS{id}{i_dat} = EEG.dipfit;
    end
    dipfitnames{id} = fitname;
end
    
save(filename,'DIPFITS','dipfitnames');
