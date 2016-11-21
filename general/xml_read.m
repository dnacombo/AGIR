function [tree, RootName, DOMnode] = xml_read(xmlfile, Pref)
%XML_READ reads xml files and converts them into Matlab's struct tree.
% 
% DESCRIPTION
% tree = xml_read(xmlfile) reads 'xmlfile' into data structure 'tree'
%
% tree = xml_read(xmlfile, Pref) reads 'xmlfile' into data structure 'tree'
% according to your preferences
%
% [tree, RootName, DOMnode] = xml_read(xmlfile) get additional information
% about XML file
%
% INPUT:
%  xmlfile	URL or filename of xml file to read
%  Pref     Preferences:
%    Pref.ItemName - default 'item' - name of a special tag used to itemize
%                    cell arrays
%    Pref.ReadAttr - default true - allow reading attributes
%    Pref.Str2Num - default true - convert strings that look like numbers 
%                   to numbers
%    Pref.NoCells - default true - force output to have no cell arrays
%
% OUTPUT:
%  tree         tree of structs and/or cell arrays corresponding to xml file
%  RootName     XML tag name used for root (top level) node
%  DOMnode      output of xmlread
%
% DETAILS:
% Function xml_read first calls MATLAB's xmlread function and than 
% converts its output ('Document Object Model' tree of Java objects) 
% to tree of MATLAB struct's. The output is often in format of nested 
% structs and cells. In the output data structure field names are based on
% XML tags, except in cases when tags produce illegal variable names.
% 
% Several special xml node types result in special tags for fields of 
% 'tree' nodes:
%  - node.CONTENT - stores data section of the node if other fields are
%    present. Usually data section is stored directly in 'node'.
%  - node.ATTRIBUTE.name - stores node's attribute called 'name'.
%  - node.COMMENT - stores node's comment section (string) .
%  - node.CDATA_SECTION - stores node's CDATA section (string).
%  - other special node types like: document fragment nodes, document type
%   nodes, entity nodes, notation nodes and processing instruction nodes
%   will be treated like regular nodes
%
% EXAMPLES:
%   See xml_examples.m
%
% See also:
%   xml_write, xmlread, xmlwrite
%
% Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com
% References:
%  - Function inspired by Example 3 found in xmlread function.
%  - Output data structures inspired by xml_toolbox structures.
  
%% default preferences
DPref.ItemName  = 'item'; % name of a special tag used to itemize cell arrays
DPref.ReadAttr  = true;   % allow reading attributes
DPref.Str2Num   = true;   % convert strings that look like numbers to numbers
DPref.NoCells   = true;   % force output to have no cell arrays
tree     = [];
RootName = [];

%% read user preferences
if (nargin>1)
  if (isfield(Pref, 'ItemName')), DPref.ItemName = Pref.ItemName; end
  if (isfield(Pref, 'ReadAttr')), DPref.ReadAttr = Pref.ReadAttr; end
  if (isfield(Pref, 'Str2Num' )), DPref.Str2Num  = Pref.Str2Num ; end
  if (isfield(Pref, 'NoCells' )), DPref.NoCells  = Pref.NoCells ; end
end

%% read xml file
try
  DOMnode = xmlread(xmlfile); 
catch
  error('Failed to read XML file %s.',xmlfile);
end

%% Find the Root node
RootNode = DOMnode.getFirstChild;
while (RootNode.getNodeType~=RootNode.ELEMENT_NODE)
  RootNode = RootNode.getNextSibling;
  if (isempty(RootNode)), return; end
end

%% parse xml file
try
  [tree RootName] = DOMnode2struct(RootNode, DPref);
catch
  error('Unable to parse XML file %s.',xmlfile);
end

  
  
%% =======================================================================
function [s sname LeafNode] = DOMnode2struct(node, Pref)
[sname LeafNode] = NodeName(node);
s = [];
  
%% === read in node data =================================================
if (LeafNode)
  if (LeafNode~=3) % supported leaf node types
    s = strtrim(char(node.getData));
    if (LeafNode==1 && Pref.Str2Num), s=str2variable(s); end
  end
  return
end
  
%% === read in children nodes ============================================
if (node.hasChildNodes)        % children present
  Child  = node.getChildNodes; % create array of children nodes
  nChild = Child.getLength;    % number of children
        
  % --- pass 1: how many children with each name ----------------------- 
  f = [];
  for iChild = 1:nChild        % read in each child
    [cname cLeaf] = NodeName(Child.item(iChild-1));
    if (cLeaf==3), continue; end % unsupported leaf node types
    if (~isfield(f,cname)), 
      f.(cname)=0;           % initialize first time I see this name
    end  
    f.(cname) = f.(cname)+1; % add to the counter
  end                        % end for iChild
  % text_nodes become CONTENT & for some reason current xmlread 'creates' a
  % lot of empty text fields so f.CONTENT value should not be trusted
  if (isfield(f,'CONTENT') && f.CONTENT>2), f.CONTENT=2; end
    
  % --- pass 2: store all the children ---------------------------------
  for iChild = 1:nChild        % read in each child
    [c cname cLeaf] = DOMnode2struct(Child.item(iChild-1), Pref);
    if (cLeaf && isempty(c))   % if empty leaf node than skip
      continue;                % usually empty text node or one of unhandled node types
    elseif (nChild==1 && cLeaf==1) 
      s=c;                     % shortcut for a common case
    else                       % if normal node
      n = f.(cname);           % how many of them in the array so far?        
      if (~isfield(s,cname))   % encountered this name for the first time
        if (n==1)              % if there will be only one of them ...
          s.(cname) = c;       % than save it in format it came in
        else                   % if there will be many of them ...
          s.(cname) = cell(1,n);
          s.(cname){1} = c;    % than save as cell array
        end
        f.(cname) = 1;         % reset the counter
      else                     % already have seen this name
        s.(cname){n+1} = c;    % add to the array
        f.(cname) = n+1;       % add to the array counter
      end  
    end
  end   % for iChild
end % end if (node.hasChildNodes)

%% === Post-processing of struct's =======================================
if (isstruct(s))
  fields = fieldnames(s);
  nField = length(fields);
    
  % --- Post-processing: convert 'struct of arrays' to 'array of struct'
  vec = zeros(size(fields)); 
  for i=1:nField, vec(i) = f.(fields{i}); end
  if (numel(vec)>1 && vec(1)>1 && var(vec)==0)    % convert from struct of
    s = cell2struct(struct2cell(s), fields, 1); % arrays to array of struct    
  end % if anyone knows better way to do above conversion please let me know.
    
  % --- Post-processing: remove special 'item' tags ---------------------
  if (isfield(s,Pref.ItemName))
    if (nField==1)
      s = s.(Pref.ItemName);         % only child: remove a level
    else
      s.CONTENT = s.(Pref.ItemName); % other children/attributes present use CONTENT
      s = rmfield(s,Pref.ItemName);
    end
  end
    
  % --- Post-processing: clean up CONTENT tags ---------------------
  if (isfield(s,'CONTENT'))
    if (iscell(s.CONTENT)) % && all(cellfun('isempty', s.CONTENT(2:end))))
      %msk = ~cellfun('isempty', s.CONTENT)
      %s.CONTENT = s.CONTENT(msk); % delete empty cells
      x = s.CONTENT;
      for i=length(x):-1:1, if ~isempty(x{i}), break; end; end
      if (i==1) 
        s.CONTENT = x{1};   % delete cell structure
      else
        s.CONTENT = x(1:i); % delete empty cells
      end
    end
    if (nField==1)
      s = s.CONTENT;      % only child: remove a level
    end
  end
end
  
%% === read in attributes ===============================================
if (node.hasAttributes && Pref.ReadAttr)
  if (~isstruct(s)),               % make into struct if is not already
    ss.CONTENT=s; 
    s=ss; 
  end  
  Attr  = node.getAttributes;     % list of all attributes
  for iAttr = 1:Attr.getLength    % for each attribute
    name  = char(Attr.item(iAttr-1).getName);  % attribute name 
    name  = genvarname_(name);    % fix name if needed
    value = char(Attr.item(iAttr-1).getValue); % attribute value
    value = str2variable(value);  % convert to number if possible
    s.ATTRIBUTE.(name) = value;   % save again
  end                             % end iAttr loop
end % done with attributes

%% === Post-processing: convert 'cells of structs' to 'arrays of structs' 
if (isstruct(s))  
  fields = fieldnames(s);     % get field names
  for iItem=1:length(s)       % for each struct in the array - usually one
    for iField=1:length(fields)
      field = fields{iField}; % get field name
      x = s(iItem).(field);   
      if (iscell(x) && all(cellfun(@isstruct,x))) % it's cells of structs  
        try                           % this operation fails sometimes 
          s(iItem).(field) = [x{:}];  % converted to arrays of structs
        catch
          if (Pref.NoCells)
            s(iItem).(field) = forceCell2Struct(x);
          end
        end % end catch
      end
    end
  end
end

%% =======================================================================
function s = forceCell2Struct(x) 
% convert cell array of structs, where not all of structs have the same
% fields, to a single array of structs
  AllFields = fieldnames(x{1});     % get field names
  CellMat = cell(length(x), length(AllFields));
  for iItem=1:length(x)      
    fields = fieldnames(x{iItem});  % get field names
    for iField=1:length(fields)
      field = fields{iField};       % get field name
      col = find(strcmp(field,AllFields),1);
      if isempty(col)
        AllFields = [AllFields; field];
        col = length(AllFields);
      end
      CellMat{iItem,col} = x{iItem}.(field);   
    end
  end
  s = cell2struct(CellMat, AllFields, 2);

%% =======================================================================
function val=str2variable(str)
% Can this string be converted to a number? if so than do it.
val = str; 
if (numel(str)==0), return; end
s = regexprep(str, '[Inf,NaN,pi,\t,\n,\d,\+,\-,\*,\.,e,i, ,E,I,\[,\],\;,\,]', '');
if (~all(~isempty(s))) 
  str(strcmp(str,'\n')) = ';'; % for parsing data tables into 2D arrays
  %str = regexprep(str, '\n', ';'); 
  num = str2num(str); % this is something created by mat2str
  if(isnumeric(num) && numel(num)>0), val=num; end
end

%% =======================================================================
function [Name LeafNode] = NodeName(node)
% get node name and make sure it is a valid variable name in Matlab.
% also get node type: 
%   LeafNode=0 - normal element node, 
%   LeafNode=1 - text node
%   LeafNode=2 -   supported non-text leaf node, 
%   LeafNode=3 - unsupported non-text leaf node
switch (node.getNodeType) 
  case node.ELEMENT_NODE
    Name = char(node.getNodeName);% capture name of the node
    Name = genvarname_(Name);     % if Name is not a good variable name - fix it  
    LeafNode = 0;
  case node.TEXT_NODE
    Name = 'CONTENT';
    LeafNode = 1;
  case node.COMMENT_NODE
    Name = 'COMMENT';
    LeafNode = 2;
  case node.CDATA_SECTION_NODE
    Name = 'CDATA_SECTION';
    LeafNode = 2;
  otherwise
    NodeType = {'ELEMENT','ATTRIBUTE','TEXT','CDATA_SECTION', ...
        'ENTITY_REFERENCE', 'ENTITY', 'PROCESSING_INSTRUCTION', 'COMMENT',...
        'DOCUMENT', 'DOCUMENT_TYPE', 'DOCUMENT_FRAGMENT', 'NOTATION'};
    Name = char(node.getNodeName);% capture name of the node
    warning('Unknown node type encountered: %s_NODE (%s)', NodeType{node.getNodeType}, Name);
    LeafNode = 3;
end

%% =======================================================================
function s = genvarname_(s)
if (length(s) > 0 && length(s) <= namelengthmax)
  if (~isempty(regexprep(s,'[A-Za-z]\w*','')))
    s = genvarname(s);  
  end
end

