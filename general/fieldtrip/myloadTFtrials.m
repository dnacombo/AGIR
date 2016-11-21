function myloadTFtrials(TFfile,Freqs)

% myloadTFtrials(TFfile,Freqs)
% load single trials TF data from TFfile.
% Freqs are the frequencies. If empty, load all frequencies. If omitted, do
% not load single trials, only average TF.

if isstruct(TFfile)
    TFfile = [TFfile.TFdatapath TFfile.TFdatafile];
end
load(TFfile,'TF');
if not(exist('Freqs','var'))
    assignin('caller','TF',TF);
    return
end
if isempty(Freqs)
    Freqs = TF.freqs;
end
Freqs = sort(Freqs);
fidx = ismember(TF.freqs, Freqs);
TF.amps = TF.amps(fidx,:,:,:);
TF.itc = TF.itc(fidx,:,:,:);
TF.freqs = Freqs;
for i_f = 1:numel(Freqs)
    load(TFfile,['TFtrials' num2str(Freqs(i_f))]);
    if i_f == 1
        ss = eval(['size(TFtrials' num2str(Freqs(i_f)) ')']);
        s = size(TF.amps,1);
        TFtrials = NaN([s ss(2:end)],'single');
    end
    eval(['TFtrials(i_f,:,:,:) = TFtrials' num2str(Freqs(i_f)) ';clear TFtrials' num2str(Freqs(i_f)) ';']);
end
assignin('caller','TFtrials',TFtrials);
assignin('caller','TF',TF);