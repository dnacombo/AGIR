function c = splitstrcell(str,cut)

if size(str,1) ~= 1
    error('str should be one line');
end
    
c{1} = str(1:cut(1));
for i = 2:numel(cut)
    if i < numel(cut)
        c{i} = str(cut(i)+1:cut(i+1));
    else
        c{i} = str(cut(i)+1:end);
    end
end