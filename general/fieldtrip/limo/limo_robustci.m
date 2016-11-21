function [mask Lo Hi]= limo_robustci(m,sd,n,bootTH1,pthresh)


if ndims(bootTH1) ~= 3
    error('Bootstrap data should have 3 dimensions')
end

nboot = size(bootTH1,3);
lo = round((nboot.*pthresh)/2);
hi = nboot - lo;

Tsorted  = sort(bootTH1,3);

Hi = m - Tsorted(:,:,hi) .* sd./sqrt(n);
Lo = m - Tsorted(:,:,lo+1) .* sd./sqrt(n);


mask = Hi >=0 | Lo <= 0;

