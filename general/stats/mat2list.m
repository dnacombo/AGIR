function list = mat2list(mat,valuename,colnames)
% list = mat2list(mat)
% 
% change a matrix to a list table with one value on each row and the
% indices of all dimensions for the corresponding element in mat.

if not(exist('valuename','var')) || isempty(valuename)
    valuename = inputname(1);
end
if not(ischar(valuename))
    error('valuename if provided should be a string')
end
letters = 'abcdefghijklmnop';
if not(exist('colnames','var')) || isempty(colnames)
    for i = 1:ndims(mat)
        colnames{i} = letters(i);
    end
end
if numel(colnames) ~= ndims(mat)
    error('column names should be same length as dimensions of mat')
end

str = '[';
for i = 1:ndims(mat)
    str = [str colnames{i} ' '];
end
str = [str '] = ind2sub(size(mat),1:numel(mat));'];
eval(str);

str = 'struct(valuename,num2cell(mat(:)),';

for i = 1:ndims(mat)
    str = [str '''' colnames{i} ''', num2cell(' colnames{i} '(:)),'];
end
str(end) = ')';
list = eval(str);

list = struct2table(list,1);


