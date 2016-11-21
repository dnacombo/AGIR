function [h him] = eeg_simplesurf(EEG,itri,h)

% eeg_simplesurf(EEG)
% display a surface image of EEG.data (averaging on the third dimension of
% EEG.data if needed).
% 
% optional inputs: eeg_simplesurf(EEG,itri,h)
% itri : trial to show (index on the third dimension of EEG.data)
% h : handle to the figure on which to plot


if nargin == 0
    EEG = evalin('caller','EEG');
end

if exist('itri','var') && ~isempty(itri)
    EEG.data = EEG.data(:,:,itri);
end
if not(exist('h','var'))
    h = figure;
    set(h,'renderer','painters');
else
    figure(h);
end
EEG.data = mean(EEG.data,3);

dat = [EEG.data];

him = imagesc([EEG.xmin:1./EEG.srate:EEG.xmax],1:EEG.nbchan,double(dat));
set(him,'tag','simplesurf')
axis xy
% view(2);
% shading flat
% axis tight
set(gca,'ytick',1:EEG.nbchan,'yticklabel', {EEG.chanlocs.labels})
c = max(abs(dat(:)));
%caxis([-c c]);
ylabel('Electrodes');
xlabel('Time (s)');
try
    title(EEG.setname);
end

