%  SAVECURRENTEDITORFILES     save the open files in the editor to matlab shortcut.
%
%  To add the current files
%  SAVECURRENTEDITORFILES                         % interactively asks you for a reference for the files
%  SAVECURRENTEDITORFILES ( 'reference' )         % adds the current files using the supplied name as reference
%  SAVECURRENTEDITORFILES ( 'reference', 1 )      % removes the files associated with the reference
%
%  The code will add a shortcut to the Matlab shortcut bar called "Saved Editor Sessions"
%
%    To retrieve the files click on the shortcut and a list will be provided where you select by entering the
%      integer associated with the files (listed numerically and using the 'reference' supplied.
%    All files are opened into the editor (All current files are closed).
%    If a file is no longer found an error message is displayed in the command window.
%
%   Note The second argument can be a structure, as shown below - this is for future customisation.
%             options.remove = 1;
%
%   Limitations:  Untested of any versions post     R2013b   (I have tried to code it so it will work by looking at the documentation)
%                 Untested on any versions prior to R2007b
%
%
%  see also edit

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Robert J Cumming: rcumming(at)matpi.com, robertjcumming(at)yahoo.co.uk
% www.matpi.com
% $Revision: 2.00 $  $Date: 25/04/2014 09:00:20 $
function SaveCurrentEditorStatus ( varargin )
%%
if nargin == 0

    SaveCurrentEditorStatus('Save',myname);
    return
end
p = fileparts(which(mfilename));
try
    load(fullfile(p,'SavedSessions.mat'));
catch
    Sessions = struct();
end
switch varargin{1}
    case 'Save'
        if nargin == 1
            myname = input ( 'Please enter a name used for future reference: ', 's' );
            if isempty ( myname )
                fprintf ( 'Warning: Nothing updated as name not provided\n' );
                return;
            end
        else
            myname = varargin{2};
        end
        currdocs = getcuropen;
        n = strcmp({Sessions.name},myname);
        if isempty(n)
            n = numel(Sessions)+1;
        end
        Sessions(n).name
        Sessions(n).docs = currdocs;
        save(fullfile(p,'SavedSessions.mat'),'Sessions');
        
    case 'Load'
        n = strcmp({Sessions.name},myname);
        if isempty(n)
            error('No session to load')
        end
        for i = 1:numel(Sessions)
            disp([num2str(i) '. ' Sessions(i).name])
        end
        myname = input ( 'Please Select (hit enter with no number to cancel):', 's' );
        if isempty ( myname )
            return;
        end
        for i = 1:numel(Sessions(n).docs)
            if exist(Sessions(n).docs(i).Filename,'file')
                edit(Sessions(n).docs(i).Filename);
            end
        end
        Docs = matlab.desktop.editor.getAll;
        for i = 1:numel(Docs)
            error set line position
        end
        
        
end
function openDocuments = getcuropen
        
switch version ( '-release' )
    case { '2011a' '2011b' '2012a' '2012b' '2013a' '2013b' '2014a' }
        openDocuments = matlab.desktop.editor.getAll;
        editorSummary = {openDocuments.Filename};
    case { '2010a' '2010b' }
        error 'not tested'
        allDocs = editorservices.getAll;
        editorSummary={allDocs.Filename};
    case { '2009a' '2009b' '2008a' '2008b' '2007a' '2007b' '2006a' '2006b' }
        error 'not tested'
        edhandle = com.mathworks.mlservices.MLEditorServices;
        editorSummaryText = char(edhandle.builtinGetOpenDocumentNames);
        sE = size(editorSummaryText);
        editorSummary = cell(sE(1),1);
        for i=1:sE(1)
            editorSummary{i} = strtrim(editorSummaryText(i,:));
        end
    otherwise
        fprintf ( 2, 'version not supported\n' );
        return
end
%%




