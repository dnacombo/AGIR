list_IFI = [];
this_block = 1;
for i = 1 : length(C{this_block}) - 1
    list_IFI = [list_IFI, C{this_block}(i+1, 6) - C{this_block}(i, 6) ]; 
end
nb_HB_skipped = find(list_IFI>1.33); %% considering if the interval between 2 flips exceeds 1.3 s, then at least one HB was skipped 
