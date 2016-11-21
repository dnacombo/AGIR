
N1 = 200;
N2 = 200;
rM =  [.85];
rP =  [.81]; 
rois = {'Confidence vs Consensus : gp 3 & 4'};
disp('==================================================================')
for i = 1:length(rM)

    r1 = []; r2 = []; z1 = []; z2 = []; Z = []; p=[];

    r1 = rM(i);
    r2 = rP(i);

    z1 = .5*log([1+r1]/[1-r1]);

    z2 = .5*log([1+r2]/[1-r2]);
    
    Z = (z1 - z2)/sqrt(1/(N1-3) +1/(N2-3));
    p = 2*(1-normcdf([abs(Z)]));

    disp([rois{i} ':   ' 'z1 = ' num2str(z1,2) ' | z2 = ' num2str(z2,2)  ...
        ' >> Z = ' num2str(Z,3) ', p = ' num2str(p,4) ])
    zvals(i) = Z;
    pvals(i) = p;
end
