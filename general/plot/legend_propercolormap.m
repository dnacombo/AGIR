function legend_propercolormap

[dum hl] = legend;
cm = colormap;
j = 1;
for i = 1:numel(hl)
    if strcmp(get(hl(i),'type'),'text')
        continue
    else
        set(get(get(hl(i),'children'),'children'),'color',cm(j,:));
        j = j+1;
    end
end

