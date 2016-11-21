function varargout = mynewlimo_resultsgui(varargin)
% MYNEWLIMO_RESULTSGUI MATLAB code for mynewlimo_resultsgui.fig
%      MYNEWLIMO_RESULTSGUI, by itself, creates a new MYNEWLIMO_RESULTSGUI or raises the existing
%      singleton*.
%
%      H = MYNEWLIMO_RESULTSGUI returns the handle to a new MYNEWLIMO_RESULTSGUI or the handle to
%      the existing singleton*.
%
%      MYNEWLIMO_RESULTSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYNEWLIMO_RESULTSGUI.M with the given input arguments.
%
%      MYNEWLIMO_RESULTSGUI('Property','Value',...) creates a new MYNEWLIMO_RESULTSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mynewlimo_resultsgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mynewlimo_resultsgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mynewlimo_resultsgui

% Last Modified by GUIDE v2.5 10-Feb-2012 17:51:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mynewlimo_resultsgui_OpeningFcn, ...
    'gui_OutputFcn',  @mynewlimo_resultsgui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mynewlimo_resultsgui is made visible.
function mynewlimo_resultsgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mynewlimo_resultsgui (see VARARGIN)

% Choose default command line output for mynewlimo_resultsgui
handles.output = hObject;

set(0,'DefaultFigureRenderer','OpenGL'); % to fix bug when setting transparency on multiple monitor config

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mynewlimo_resultsgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mynewlimo_resultsgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushload.
function pushload_Callback(hObject, eventdata, handles)
% hObject    handle to pushload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

EEG = pop_loadset;
if isempty(EEG)
    try
        EEG = handles.EEG;
    catch
        return
    end
end
EEG.times = EEG.xmin:1/EEG.srate:EEG.xmax;
EEG.data = mean(EEG.data,3);
handles.EEG = EEG;
try
    cd(EEG.filepath);
    set(handles.edit_dir,'String',EEG.filepath);
end
txt = {};
if not(isempty(EEG))
    txt{end+1} = ['Setname: ' EEG.setname];
    txt{end+1} = ['File: ' EEG.filename];
    txt{end+1} = [num2str(EEG.nbchan) 'x' num2str(EEG.trials) 'x' num2str(EEG.pnts) ' data'];
    txt{end+1} = [num2str(EEG.xmin) '-' num2str(EEG.xmax) ' s (' num2str(EEG.srate) ' Hz)'];
    if EEG.trials>1
        txt{end+1} = '';
        txt{end+1} = ['Warning : we have several trials here. They are averaged.'];
        txt{end+1} = '';
    end
    txt{end+1} = '_______________________________________';
    txt{end+1} = EEG.history;
end
set(handles.text_loaded,'String',txt);

guidata(hObject,handles);

push_eeg_simplesurf_Callback(hObject, eventdata, handles)


% --- Executes on button press in push_pop_topoplot.
function push_pop_topoplot_Callback(hObject, eventdata, handles)
% hObject    handle to push_pop_topoplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pop_topoplot(handles.EEG);

% --- Executes during object creation, after setting all properties.
function push_pop_topoplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to push_pop_topoplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in push_eeg_simplesurf.
function push_eeg_simplesurf_Callback(hObject, eventdata, handles)
% hObject    handle to push_eeg_simplesurf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hsimplesurffig = figure(63854);
set(handles.hsimplesurffig,'numbertitle','off');

handles.currentelec = 1;

% topo
handles.htopofig = figure(64822);clf;
set(handles.htopofig,'numbertitle','off', 'defaultaxesfontsize',14);
topoplot(handles.EEG.data(:,handles.currentelec),handles.EEG.chanlocs,'conv','on','emarker2',{handles.currentelec,'.','k',26,1});
title({handles.EEG.chanlocs(1).labels [num2str(handles.EEG.xmin,'%.3g') 's']})
colorbar

% time course
handles.hplotfig = figure(6416584);clf;
set(handles.hplotfig,'numbertitle','off', 'defaultaxesfontsize',14);
set(gca,'ylimmode','auto');
handles.hplot = plot(handles.EEG.times,handles.EEG.data(handles.currentelec,:));
set(handles.hplot,'Tag','timecourse');
yl = ylim;xlim([handles.EEG.xmin handles.EEG.xmax]);
hold on
plot([handles.EEG.times(1) handles.EEG.times(1)],yl,':k');
hold off
title(handles.EEG.chanlocs(handles.currentelec).labels)

% simple surf
figure(handles.hsimplesurffig);
[handles.hsimplesurffig handles.hsimplesurf] = eeg_simplesurf(handles.EEG,[],handles.hsimplesurffig);
set(handles.hsimplesurf,'Tag','simplesurf');
set(gcf,'UserData',[handles.hplotfig handles.htopofig],'closerequestfcn','try close(get(gcbo,''UserData''));end; delete(gcbo)');
if all(handles.EEG.data(:) >= 0)
    colormap hot
else
    colormap jet
end
cx = caxis;
caxis([-max(abs(cx)) max(abs(cx))]);
colorbar
figure(handles.hsimplesurffig);
handles.datatipsimplesurf = datacursormode(handles.hsimplesurffig);
set(handles.datatipsimplesurf,'UpdateFcn',{@updatERP, handles.figure1});
datacursormode on
guidata(hObject,handles);

function output_txt = updatERP(obj,eventobj,hmainfig)

handles = guidata(hmainfig);

EEG = handles.EEG;
pos = get(eventobj,'Position');

[handles.currentelec] = chnb(pos(2));

[handles.currentelec elec] = chnb(handles.currentelec);

output_txt = {['time: ',num2str(pos(1),4)],...
     ['elec: ',elec{1}]};
figure(handles.hplotfig)
clf
set(gcf,'name',['Time course at ' elec{1}])
plot(EEG.times,EEG.data(handles.currentelec,:));
if get(handles.check_robust,'Value')
    updatestatview(handles);
else
    set(handles.hsimplesurf,'AlphaData',1)
end
yl = ylim;xlim([EEG.xmin EEG.xmax]);
hold on
plot([pos(1) pos(1)],yl,':k');
hold off
title(elec{1})

figure(handles.htopofig)
clf
set(gcf,'name',['Topographical map at ' num2str(pos(1),4) 's'] )
topoplot(EEG.data(:,EEG.times == pos(1)),EEG.chanlocs,'conv','on','emarker2',{handles.currentelec,'.','k',26,1});
title({elec{1} [num2str(pos(1),4) 's']})
colorbar
figure(handles.hsimplesurffig);
guidata(handles.figure1,handles);

function edit_dir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dir as text
%        str2double(get(hObject,'String')) returns contents of edit_dir as a double


% --- Executes during object creation, after setting all properties.
function edit_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',cd)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_dir.
function edit_dir_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

p = uigetdir(cd,'Choose a directory');
try
    cd(p);
    set(hObject,'String',p)
end


% --- Executes on button press in push_robustCI.
function push_robustCI_Callback(hObject, eventdata, handles)
% hObject    handle to push_robustCI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    fprintf('Loading saved confidence intervals... ')
    s = warning('off','MATLAB:load:variableNotFound');
    load(handles.EEG.statsfile,'CIpthresh');
    warning(s.state,s.identifier);
    if CIpthresh == str2num(get(handles.edit_pthresh,'String'))
        handles.stats.CIpthresh = CIpthresh;
        s = warning('off','MATLAB:load:variableNotFound');
        load(handles.EEG.statsfile,'CImask','CILo','CIHi');
        warning(s.state,s.identifier);
        handles.stats.CImask = CImask;
        handles.stats.CILo = CILo;
        handles.stats.CIHi = CIHi;
        fprintf('done\n')
    else
        handles.stats = struct;
        error
    end
catch
    fprintf('failed\n')
    fprintf('Computing robust confidence intervals... ')
    handles.stats.CIpthresh = str2num(get(handles.edit_pthresh,'String'));
    stats = load(handles.EEG.statsfile,'m','sd','n','boottH0');
    [handles.stats.CImask handles.stats.CILo handles.stats.CIHi]  = ...
        limo_robustci(stats.m,stats.sd,stats.n,stats.boottH0,handles.stats.CIpthresh);
    struct2ws(handles.stats);
    save(handles.EEG.statsfile,'-append','CIpthresh','CImask','CILo','CIHi');
    fprintf('done\n')
end

guidata(hObject,handles);

% --- Executes on button press in push_robustSTAT.
function push_robustSTAT_Callback(hObject, eventdata, handles)
% hObject    handle to push_robustSTAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    fprintf('Loading saved robust cluster mask... ')
    s = warning('off','MATLAB:load:variableNotFound');
    load(handles.EEG.statsfile,'STATpthresh');
    warning(s.state,s.identifier);
    if STATpthresh == str2num(get(handles.edit_pthresh,'String'))
        handles.stats.STATpthresh = STATpthresh;
        s = warning('off','MATLAB:load:variableNotFound');
        load(handles.EEG.statsfile,'STATmask');
        warning(s.state,s.identifier);
        handles.stats.STATmask = STATmask;
        fprintf('done\n')
    else
        handles.stats = struct;
        error
    end
catch
    fprintf('failed\n')
    handles.stats.STATpthresh = str2num(get(handles.edit_pthresh,'String'));
    fprintf('Loading Bootstrap data... ')
    stats = load(handles.EEG.statsfile,'t','p','n','boottH0','bootpH0');
    fprintf('done\n')
    try
        fprintf('Loading channel neighborhood matrix... ')
        s = warning('off','MATLAB:load:variableNotFound');
        load(handles.EEG.channeighbstructmatfile)
        warning(s.state,s.identifier);
        handles.stats.channeighbstructmat = channeighbstructmat;
        fprintf('done\n')
    catch
        fprintf('failed.\n')
        resp = questdlg('Create a new channel neighborhood matrix?','No channel neighborhood file name in EEG structure','Yes','No, load one','Cancel','Cancel');
        switch resp
            case 'Yes'
                limo_expected_chanlocs
                load expected_chanlocs
                handles.stats.channeighbstructmat = channeighbstructmat;
            case 'No, load one'
                [f p] = uigetfile('*.mat','Choose an expected channel location file.','expected_chanlocs.mat');
                load([p f]);
                channeighbstructmatfile = [p f];
                handles.EEG.channeighbstructmatfile = channeighbstructmatfile;
                handles.stats.channeighbstructmat = channeighbstructmat;
                pop_saveset(handles.EEG,'savemode','resave');
            case 'Cancel'
                return
        end
    end
    fprintf('Computing cluster size correction for multiple comparisons... ')
    handles.stats.STATmask = limo_clusterthresh(stats.t,stats.p,stats.boottH0, stats.bootpH0,handles.stats.STATpthresh,handles.stats.channeighbstructmat);
    struct2ws(handles.stats);
    save(handles.EEG.statsfile,'-append','STATpthresh','STATmask');
    fprintf('done\n')
end

guidata(hObject,handles);

function edit_pthresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pthresh as text
%        str2double(get(hObject,'String')) returns contents of edit_pthresh as a double
set(hObject,'String',sprintf('%0.2g',str2num(get(hObject,'String'))));

% --- Executes during object creation, after setting all properties.
function edit_pthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function struct2ws(s,varargin)

% struct2ws(s,varargin)
%
% Description : This function returns fields of scalar structure s in the
% current workspace
% __________________________________
% Inputs :
%   s (scalar structure array) :    a structure that you want to throw in
%                                   your current workspace.
%   re (string optional) :          a regular expression. Only fields
%                                   matching re will be returned
% Outputs :
%   No output : variables are thrown directly in the caller workspace.
%
% _____________________________________
% See also : ws2struct ; regexp
%
% Maximilien Chaumon v1.0 02/2007


if length(s) > 1
    error('Structure should be scalar.');
end
if not(isempty(varargin))
    re = varargin{1};
else
    re = '.*';
end

vars = fieldnames(s);
vmatch = regexp(vars,re);
varsmatch = [];
for i = 1:length(vmatch)
    if isempty(vmatch{i})
        continue
    end
    varsmatch(end+1) = i;
end
for i = varsmatch
    assignin('caller',vars{i},s.(vars{i}));
end


% --- Executes on button press in check_robust.
function check_robust_Callback(hObject, eventdata, handles)
% hObject    handle to check_robust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_robust

handles.transp = .3;
guidata(hObject,handles);

if isfield(handles,'hsimplesurf') && ishandle(handles.hsimplesurf)
    updatestatview(handles)
end


function updatestatview(handles)

set(get(handles.hsimplesurf,'parent'),'color','k')
set(handles.hsimplesurf,'AlphaData',handles.stats.STATmask*handles.transp+1-handles.transp)
figure(handles.hplotfig);
hold on;
try delete(handles.hplotarea); end
handles.hplotarea = fill([handles.EEG.times handles.EEG.times(end:-1:1)], ...
    [handles.stats.CILo(handles.currentelec,:) handles.stats.CIHi(handles.currentelec,end:-1:1)], ...
    'r', 'EdgeColor', 'none', 'FaceColor', [1 0 0 ], 'facealpha', .1);
hold off;
