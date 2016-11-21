function [blocks_def] = identifying_blocks(list_trig_for_blocks)

j=1;
blocks_def = [];
for i = 1 : 2: length(list_trig_for_blocks)
    a = list_trig_for_blocks(i);
    b = list_trig_for_blocks(i+1);
    blocks_def(j,:) = [a b];    
        j=j+1;
end

