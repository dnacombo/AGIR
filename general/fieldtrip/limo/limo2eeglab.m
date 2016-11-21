function EEG = limo2eeglab(limofile)

% EEG = limo2eeglab(limofile)
% 
% Put a data matrix contained in limofile (2-3D) in place of the data in an EEG structure/file.
% update the structure fields minimally.
%

EEG = eeg_emptyset;
if isstr(limofile)
    p = fileparts(limofile);
    if isempty(p)
        p = cd;
    end
    tmp = load(limofile);
    f = fieldnames(tmp);
    limodat = tmp.(f{1});
    clear tmp
else
    limodat = limofile;
    limofile = sprintf('<[%g %g %g] data>',size(mat));
    p = cd;
end
try
    load(fullfile(p,'LIMO.mat'));
catch
    error('No LIMO.mat found. I need it.')
end
if not(any(ndims(limodat) == [2 3]))
    error('Limo data should have 2 or 3 dimensions')
end
EEG.data = limodat;
EEG.chanlocs = LIMO.data.chanlocs;
EEG.srate = LIMO.data.sampling_rate;
EEG.xmin = LIMO.data.start + (LIMO.data.trim1-1)/LIMO.data.sampling_rate;
EEG.setname = sprintf('%s',limofile);
EEG.comments = sprintf('limofile: %s',fullfile(p,limofile));
EEG.history = [sprintf('EEG = limo2eeglab(''%s'')',limofile)];
EEG.filename = limofile;
EEG.datfile = limofile;

EEG = eeg_checkset(EEG);
