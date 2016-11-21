function [Im_pwr2s,Im_fft2_ls] = compute_2dpwr(Im,showresult,logtrans)

% [Im_pwr2s,Im_fft2_ls] = compute_2dpwr(Im,showresult,logtrans)
% 




fs = size(Im,1);

PQ = paddedsize(size(Im),'pwr2');
Im_fft2 = fft2(double(Im));

Im_fft2_ls= fftshift(log(1+abs( ( fft2( double(Im),fs,fs) ) )));

Im_pwr2 = (conj(Im_fft2)/ fs).*Im_fft2;
Im_pwr2_log = log(1+abs(Im_pwr2));
imsize = size(Im_pwr2);

if logtrans

    Im_pwr2s = [Im_pwr2_log(imsize(1)/2+1:imsize,1:imsize(1)/2); ...
        Im_pwr2_log(1:imsize(1)/2,1:imsize(1)/2)];
else
    %       Im_pwr2s = [Im_pwr2(imsize(1)/2+1:imsize,1:imsize(1)/2); ...
    %           Im_pwr2(1:imsize(1)/2,1:imsize(1)/2)];

    Im_pwr2s = [Im_pwr2(imsize(1)/2+1:imsize, 1:imsize(1)/2); ...
        Im_pwr2(1:imsize(1)/2,1:imsize(1)/2)];
    % 1:63, 1:63
end

if showresult
    figure(18);clf;
    set(gcf,'colormap',gray,...
        'units','pixels','pos',[300 200 256 256]);
    imagesc(Im_fft2_ls)
    axis square

    figure(19),
    clf;
    if logtrans
        imagesc(Im_pwr2s,[0 30]);
    else
        imagesc(Im_pwr2s,[0 1e6]);
    end
    set(gcf,'colormap',gray,...
        'units','pixels',...
        'pos',[360 400 256 256],...
        'colormap',jet);
    colorbar
    axis square
end
