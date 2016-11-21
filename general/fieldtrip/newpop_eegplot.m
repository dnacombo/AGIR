function newpop_eegplot(EEG, icacomp, superpose, reject ,varargin)

% same as pop_eegplot except that you can select single electrodes on
% individual trials.

addpath('/home/chaumonm/Dropbox/MATLAB/general/eeglab13');
% because I had to correct a bug in eegplot2trial. Bug #1179

if nargin < 1
	help pop_eegplot;
	return;
end;	
if nargin < 2
	icacomp = 1;
end;	
if nargin < 3
	superpose = 0;
end;
if nargin < 4
	reject = 0;
end;

myoptions = {'ctrlselectcommand'
    {'eegplot_selectelec(gcbf);'  'eegplot(''defmotioncom'', gcbf);'    'eegplot(''defupcom'',     gcbf);'};
    };
topcommand = 'doit = getpref(''warnnewpopeegplot'',''doit'',1);if doit;uiwait(warndlg(''As of 18.03.15, this function has changed. It does NOT reject and interpolate trials anymore. Instead, you need to do it yourself in your script, by entering e.g. the text shown in the command window now. This message will show up only once.''));end;fprintf(''############### Text to copy ##############\nkeyboard;\n[EEG, com] = pop_selectiveinterp(EEG);\nEEG = eegh(com,EEG);\n[EEG, com] = pop_rejepoch(EEG, find(EEG.reject.rejmanual), 1);\nEEG = eegh(com,EEG);\n################### end text to copy ##############\n'');setpref(''warnnewpopeegplot'',''doit'',0);';

eeglab redraw
pop_eegplot(EEG, icacomp, superpose, reject ,topcommand,varargin{:},myoptions{:});





