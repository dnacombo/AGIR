function mask = limo_clusterthresh(T,P,bootTH0,bootPH0,pthresh, neighbors)
% This implements a cluster size correction for multiple comparisons
% mask = clusterthresh(T,P,bootTH0,bootPH0,pthresh, neighbors)

nboot = size(bootTH0,3);
U = round((1-pthresh)*size(bootTH0,ndims(bootTH0)));

if size(bootTH0,1)>1 % many electrodes
    minnbchan = 2; % we take at leat two channels to make a cluster
    % compute bootstrap clusters
    boot_maxclustersum=zeros(nboot,1); 
    for s=1:nboot
        boot_maxclustersum(s) = limo_getclustersum(bootTH0(:,:,s).^2,bootPH0(:,:,s),neighbors,minnbchan,pthresh);
    end
    sort_boot_maxclustersum = sort(boot_maxclustersum,1);
    mask = limo_cluster_test(T.^2,P,sort_boot_maxclustersum(U),neighbors,minnbchan,pthresh);
elseif size(bootTH0,1)==1 % one electrode
    th = limo_ecluster_make( squeeze(bootTH0).^2,squeeze(bootPH0),pthresh );
    sigcluster = limo_ecluster_test( T.^2,P,th,pthresh );
    mask = sigcluster.elec;
end


