function xlswrite2(ms,sheader,scolnames,filename,sheets)
% xlswrite2     Easily create SEVERAL Excel spreadsheets from MATLAB
%
% Adapted from xlswrite
%
%
%  xlswrite(m,header,colnames,filename,sheets) creates a Microsoft Excel workbook using
%  the MATLAB ActiveX interface.  Microsoft Excel is required.
%
%Inputs:
%    m          Cellarray of Matrices to write to file
% (Optional):
%    header     Cellarray of header information. Use cell arrays of cell arrays for multiple lines
%                  DO NOT USE multiple row character arrays!!
%    colnames   Cellarray of (Cell array of strings) Column headers.  One cell element per column
%    filename   (string) Name of Excel file.  If not specified, contents will
%                  be opened in Excel.
%    sheetname:  Cellarray of Names of sheets to write data to . The
%                defaults are 'Feuil1' 'Feuil2' 'Feuil3' ...
%                   Non existing sheets are created.
%

% Scott Hirsch
% The MathWorks
% This is provided free, no warranty, ...
%
% Copied from ActiveX example in documentation
%
% In collaboration with:
% Fahad Al Mahmood fahad@al-mahmood.com
% Fahad developed the capabilities for writing to multiple sheets
%
% Dragon.  dragon5645995@sina.com.cn
% Dragon fixed a bug when writing out 52 or more columns
% 
% Maximilien Chaumon.       chaumon@chups.jussieu.fr
% extended xlswrite to write several spreadsheets and
% create it if needed.

% Parse inputs
nsheets = length(ms);
if nargin<2 | isempty(sheader)
    sheader = cell(1,nsheets);
end;

if nargin<3 | isempty(scolnames)
    scolnames = cell(1,nsheets);
end;

if nargin<4 | isempty(filename)
    visible = 1;    % Not saving to a file, so make Excel visible
    filename = '';
else
    visible = 0;    % Saving to a file.  Keep Excel hidden
end;

if nargin < 5 | isempty(sheets)
    for i = 1:nsheets
        sheets{i} = ['Feuil' num2str(i)];
    end
end;


for i = 1:nsheets
    [nr,nc] = size(ms{i});
    if nc>256
        error(['Sheet ' sheets{i} ': Matrix is too large.  Excel only supports 256 columns']);
    end;
end

% check if extension in filename
[dum,fn,ex] = fileparts(filename);
if isempty(ex)
    filename = [filename '.xls'];
end
% try
% Open Excel, add workbook, change active worksheet, 
% get/put array, save.
% First, open an Excel Server.
Excel = actxserver('Excel.Application');

% Three cases:
% * Open a new workbook, but don't save (filename is empty)
% * Open a new workbook, save with given file name
% * Open an existing workbook
if isempty(filename)
    % Insert a new workbook.
    op = invoke(Excel.Workbooks, 'Add');
    
elseif exist(filename,'file')==0
    % The following case if file does not exist (Creating New File)
    op = invoke(Excel.Workbooks,'Add');
    if regexp(filename,'\<\w:.*')
        invoke(op, 'SaveAs', filename);
    else
        invoke(op, 'SaveAs', [pwd filesep filename]);
    end
    new=1;
else
    % The following case if file does exist (Opening File)
    disp(['Opening Excel File ...(' filename ')']);
    if regexp(filename,'\<\w:.*')
        op = invoke(Excel.Workbooks, 'open', filename);
    else
        op = invoke(Excel.Workbooks, 'open', [pwd filesep filename]);
    end
    new=0;
end

%If the user does not specify a filename, we'll make Excel visible
%If they do, we'll just save the file and quit Excel without ever making
% it visible
set(Excel, 'Visible', visible);      %You might want to hide this if you autosave the file


for i_sheet = 1:nsheets
    m = ms{i_sheet};
    sheetname = sheets{i_sheet};
    colnames = scolnames{i_sheet};
    header = sheader{i_sheet};
    % Make the specified sheet active.
    try 
        Sheets = Excel.ActiveWorkBook.Sheets;
        target_sheet = get(Sheets, 'Item', sheetname);
    catch
        % Error if the sheet doesn't exist. 
        % Create it.
        % The alternative to try/catch is to call xlsfinfo to see if the sheet exists, but
        % that's really slow.
        ns = invoke(Sheets,'Add');
        set(ns,'Name',sheetname);
        target_sheet = get(Sheets, 'Item', sheetname);
    end;
    
    invoke(target_sheet, 'Activate');
    
    [nr,nc] = size(m);
    if nc>256
        error('Matrix is too large.  Excel only supports 256 columns');
    end;
    
    
    
    %Write header
    Activesheet = Excel.Activesheet;
    if isempty(header)
        nhr=0;
    elseif iscell(header)
        nhr = length(header);       %Number header rows
        for ii=1:nhr
            ActivesheetRange = get(Activesheet,'Range',['A' num2str(ii)],['A' num2str(ii)]);
            set(ActivesheetRange, 'Value', header{ii});
        end;
    else
        nhr = 1;                   %Number header rows
        ActivesheetRange = get(Activesheet,'Range','A1','A1');
        set(ActivesheetRange, 'Value', header);
    end;
    
    
    %Add column names
    
    if ~isempty(colnames)
        nhr = nhr + 1;      %One extra column name
        ncolnames = length(colnames);
        for ii=1:ncolnames
            colname = localComputLastCol('A',ii);
            %    cellname = [char(double('A')+ii-1) num2str(nhr+1)];
            cellname = [colname num2str(nhr)];
            ActivesheetRange = get(Activesheet,'Range',cellname,cellname);
            set(ActivesheetRange, 'Value', colnames{ii});
        end;
    end;
    
    
    % Put a MATLAB array into Excel.
    FirstRow = nhr+1;           %You can change the first data row here.  I start right after the headers
    if nr == 0
        LastRow = FirstRow;
    else
        LastRow = FirstRow+nr-1;
    end
    FirstCol = 'A';         %You can change the first column here
    LastCol = localComputLastCol(FirstCol,nc);
    ActivesheetRange = get(Activesheet,'Range',[FirstCol num2str(FirstRow)],[LastCol num2str(LastRow)]);
    set(ActivesheetRange, 'Value', m);
    
end % i_sheet

% If user specified a filename, save the file and quit Excel

% If user specified a filename, save the file and quit Excel
if ~isempty(filename)
    [pathstr,name,ext] = fileparts(filename);
    if isempty(pathstr)
        pathstr = pwd;
    end;
    
    invoke(op, 'Save');
    %     invoke(Workbook, 'SaveAs', [pathstr filesep name ext]);
    invoke(Excel, 'Quit');
    
    [pathstr,name,ext] = fileparts(filename);
    disp(['Excel file ' name '.xls has been created.']);
end;

%Delete the ActiveX object
delete(Excel)
% catch
%     invoke(Excel, 'Quit');
%     delete(Excel);
%     error(lasterr);
% end
function LastCol = localComputLastCol(FirstCol,nc);
% Comput the name of the last column where we will place data
%Input
%  FirstCol  (string) name of first column
%  nc        total number of columns to write

%Excel's columns are named:
% A B C ... A AA AB AC AD .... BA BB BC ...
FirstColOffset = double(FirstCol) - double('A');    %Offset from column A
if nc<=26-FirstColOffset       %Easy if single letter
    %Just convert to ASCII code, add the number of needed columns, and convert back
    %to a string
    if nc == 0
        LastCol = FirstCol;
    else
        LastCol = char(double(FirstCol)+nc-1);
    end
else
    % Fix for 52 or more columns generously provided by dragon5645995@sina.com.cn
    ng = ceil(nc/26);       %Number of groups (of 26)
    rm = rem(nc,26)+FirstColOffset;        %How many extra in this group beyond A
    if rem(nc,26)==0
        rm=26;
    end
    LastColFirstLetter = char(double('A') + ng-2);
    LastColSecondLetter = char(double('A') + rm-1);
    LastCol = [LastColFirstLetter LastColSecondLetter];
end;