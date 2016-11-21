function S = txt2struct(txt)
% S = txt2struct(txt)
% change txt table into a structure array.
% assumes row 1 of txt is valid field names
% fields are filled with each row element.

fs = txt(1,:);txt(1,:) = [];
fs = regexprep(fs,'[^a-zA-Z]','');
s = 'struct(';
for i = 1:numel(fs)
    if isempty(fs{i})
        fs{i} = ['field' num2str(i,'%02d')];
    end
    s= [s '''' fs{i} ''',[],'];
end
s(end) = ')';
S = eval(s);

for i = 1:size(txt,1)
    for j = 1:size(txt,2)
        S(i).(fs{j}) = txt{i,j};
    end
end
