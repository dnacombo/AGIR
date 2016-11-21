function liste = toksinregexp(re)
% r�cup�re les tokens et renvoit une structure int�ressante dans une
% expression r�guli�re


if nargin == 0
    re = 'sujet(01|02|03|04|05|06|07|08|09|10|11|12)_repet(2-5|9-12)_COND_(P|nP).*evoked\.minf$';
end
d = '(';
f = ')';
ou = '|';
ext = regexp(re,'\\.\w*?\$?$');
if not(isempty(ext))
    ext = ext-1;
end

pos_d = regexp(re,['\' d]);
pos_f = regexp(re,['\' f]);
pos_d_names = [1 pos_f+1];
pos_f_names = [pos_d - 1 ext];

if isempty(pos_d) liste = {}; return; end
if not(length(pos_d) == length(pos_f)) error(['unbalanced ' d ' ' f '.']); end

for i_d = 1:length(pos_d)
    if pos_d(i_d) == 1 || ( strcmp(re(pos_d(i_d)-1),'_') && ~strcmp(re(pos_f(i_d)+1),'_') )
        % after
        tiks{i_d}.name = re(pos_d_names(i_d+1) : pos_f_names(i_d+1));
    elseif pos_f(i_d) == length(re) || ( ~strcmp(re(pos_d(i_d)-1),'_') && strcmp(re(pos_f(i_d)+1),'_') )
        % before
        tiks{i_d}.name = re(pos_d_names(i_d) : pos_f_names(i_d));
    elseif strcmp(re(pos_d(i_d)-1),'_') && strcmp(re(pos_f(i_d)+1),'_')
        % no
        tiks{i_d}.name = '';
    else
        error('wtf');
    end
    
    if isempty(tiks{i_d}.name)
        tiks{i_d}.name = 'no name';
    end
    tiks{i_d}.name = strrep(tiks{i_d}.name,'_','');

    toks{i_d} = re(pos_d(i_d) + 1 : pos_f(i_d) - 1);
    pos_ou = [0 regexp(toks{i_d},['\' ou]) length(toks{i_d}) + 1];
    for i_ou = 1 : length(pos_ou) - 1
        tiks{i_d}.toks{i_ou}= toks{i_d}(pos_ou(i_ou) + 1 : pos_ou(i_ou + 1)-1);
    end
end

liste = tiks;
