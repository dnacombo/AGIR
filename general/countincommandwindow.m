function countincommandwindow(c)

persistent h t

if not(isempty(h)) && now - t < 5
    fprintf(repmat('\b',1,h));
end
h = fprintf('%g',c);
t = now;