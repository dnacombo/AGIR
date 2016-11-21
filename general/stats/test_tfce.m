% addpath('/home/chaumonm/matlab/general')

data = spatialPatternn([1,500,300],-3);

tt = tfce(data);

figure(222);
subplot(2,2,1)
h = imagesc(squeeze(data));
title('data')
subplot(2,2,2)
im = zeros(size(data));
steps = linspace(min(data(:)),max(data(:)),4);
for i = 1:numel(steps)
    im = im + double(bwlabeln(data > steps(i)) ~= 0);
end
h = imagesc(squeeze(im));
title('fixed threshold clusters')
subplot(2,2,3)
h = imagesc(squeeze(tt));
title('TFCE clusters')
subplot(2,2,4)
h = imagesc(squeeze(data));
set(h,'alphadata',squeeze(tt));
title('data masked by TFCE clusters')


%% now test tfce in clusterstats
s = [64,500,20];
null = zeros(s);
for i = 1:s(3)
    if not(rem(i,100))
        disp(i)
    end
    null(:,:,i) = spatialPatternn(s(1:2),-1);
end
%%
one = spatialPatternn([64,500],-1);
one(:,250:313) = one(:,250:313) + .05 * gabor(64);
%%
one = (one - min(one(:)))/range(one(:));
%%
figure;
subplot(2,2,1)
imagesc(one);
%%
tt = tfce(reshape(one,[1,s(1:2)]),[],.01);
subplot(2,2,3)
imagesc(squeeze(tt))
%%
[msk,pct] = clusterstats(reshape(null,[1,s]),[],reshape(one,[1,s(1:2)]),[],[],.05,[],1);
msk = reshape(msk,s(1:2));
subplot(2,2,2)
imagesc(double(squeeze(msk)))
%%
subplot(2,2,4)
h = imagesc(one);
set(h,'alphadata',.5*double(squeeze(msk))+.5)






