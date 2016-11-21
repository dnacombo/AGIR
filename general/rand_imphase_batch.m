Im=[];Img=[];swap_Im=[];ImFourier =[]; Amp=[];Phase=[];SwapPhase=[];ImSwp=[];

%  imagename = 'baseball_bat1';
%  imagename = 'newlabpic';
%   imagename = 'monitor_onwhitebg';
%   imagename = 'gslh';
%  imagename = 'PC_keybd_clean';
% basedir= './filt_uf/';
%    basedir= '/homes/11/kestas/images/cibustimuli/new/sparse/';
%       basedir= '/homes/11/kestas/images/cibustimuli/new/Dense/';
      basedir = '/space/bomba/3/users/kestas/cibustimuli/';
      
    newdirname = '/space/bomba/3/users/kestas/cibustimuli/phase_scrambled_DA_feb9_2009/';
     
    if ~exist('newdirname')
    mkdir(newdirname)
    end
    
    
  images =dir([basedir upper('*jpg')]) 
  
  for i=1:length(images)
  imagename = images(i).name;
 
  img = imread(...
     [basedir imagename]);
% Im = imread(...
% [basedir imagename '.jpg']);
% %    './pilot_isoScone/match_scenes/kitchen1.jpg'); %image to be manipulated
% 
% swap_Im = imread(...
%     [basedir imagename '.jpg']);
    %'./pilot_isoScone/match_scenes/kitchen2.jpg'); %image whose phase will be swapped in
% swap_Im = rand(size(Im,1), size(Im,2));

Im = img(:,:,1);
swap_Im = Im;

% Img = mat2gray(double(Im));
Img = double(Im);
meanint = mean(mean(Img));
swap_Im = double(swap_Im);
% swap_Im = mat2gray(double(swap_Im));
%read and rescale (0-1) image
lpcutoff = 6;
bpcutoff1 = 1;
bpcutoff2 = 60; %round(min(size(Im))/2 - 1);
hpcutoff = 40;
filtertype = 'lp';
filtername = 'gaussian';

L = lpfilter(filtername, size(swap_Im,1), size(swap_Im,2), lpcutoff,1);
H = hpfilter(filtername, size(swap_Im,1), size(swap_Im,2), hpcutoff,1);
ImSize = size(Img);




if numel(ImSize)==2
    numLayers = 1;
else
    numLayers = 3;
end

for layer = 1:numLayers
    ImFourier(:,:,layer) = fft2(Img(:,:,layer));       
%Fast-Fourier transform
    Amp(:,:,layer) = abs(ImFourier(:,:,layer));       
%amplitude spectrum
    Phase(:,:,layer) = angle(ImFourier(:,:,layer));   
%phase spectrum3
%     Phase(:,:,layer) = Phase(:,:,layer) + SwapPhase;
%add random phase to original phase
if strcmp(filtertype, 'lp')
    disp('lp-filtering')
    fswap_Im = dftfilt(swap_Im, L);
    filtPhase = angle(fft2( fswap_Im(:,:,1) ));
elseif strcmp(filtertype, 'hp')
    fswap_Im = dftfilt(swap_Im, H);
    fswap_Im = fswap_Im+meanint;
    filtPhase = angle(fft2( fswap_Im(:,:,1) ));
    disp('hp-filtering')

elseif strcmp(filtertype, 'bp')
    BP = L + H;
%     BP = BP-.5;
    fswap_Im = dftfilt(swap_Im, BP);
    fswap_Im = fswap_Im;
%     filtPhase = angle(fft2( fswap_Im(:,:,1) ));
    disp('bandpass-filtering')
else
disp('no filtering')
    filtPhase = angle(fft2( swap_Im(:,:,1) ));
end

%swap phase of Img and swap_Im
% Phase(:,:,layer) = Phase+filtPhase;
rPhase = fftshift(Phase);
sizerow = size(Img,1);
sizecol = size(Img,2);


lprow1 = round(sizerow/2)-bpcutoff2;
lprow2 = round(sizerow/2)+bpcutoff2;
lpcol1 = round(sizecol/2)-bpcutoff2;
lpcol2 = round(sizecol/2)+bpcutoff2;

sizefiltrow = (lprow2-lprow1)+1;
sizefiltcol = (lpcol2-lpcol1)+1;

irow1 = round(sizerow/2)-bpcutoff1;
irow2 = round(sizerow/2)+bpcutoff1;
icol1 = round(sizecol/2)-bpcutoff1;
icol2 = round(sizecol/2)+bpcutoff1;

iPhase = rPhase(irow1:irow2,icol1:icol2);

lpPhase = Phase(lprow1:lprow2,lpcol1:lpcol2);
rlpPhase = lpPhase;
zcross = 0;

for m = 1:10 %make 10 phasescram version of image
    
for r=1:round(sizefiltrow/1)
    randrow = randperm(sizefiltcol);
rlpPhase(r,:) = lpPhase(r,randrow);

for c = 1:round(sizefiltcol/1)

%      if  ( sign(rlpPhase(r,c)) ~= sign(lpPhase(r,c)))
        zcross = zcross+1;
         if rand < 0.9999
        if ( sign(rlpPhase(r,c)) < 0 )
            rlpPhase(r,c) = rlpPhase(r,c)-(2*pi);
        elseif (sign(rlpPhase(r,c)) > 0 )
            rlpPhase(r,c) = rlpPhase(r,c)+(2*pi);
        end
          end
%      end
%         
end

end

disp( ['Found ' num2str(zcross) ' zero crossings' ] )
% rlpPhase(irow1:irow2,icol1:icol2) = iPhase;
orow1 = round(sizefiltrow/2)-bpcutoff1;
orow2 = round(sizefiltrow/2) + bpcutoff1;

ocol1 = round(sizefiltcol/2)-bpcutoff1;
ocol2 = round(sizefiltcol/2) + bpcutoff1;


rlpPhase(orow1:orow2,ocol1:ocol2)=iPhase;

rPhase(lprow1:lprow2, lpcol1:lpcol2) = rlpPhase;
% rPhase(lprow1:lprow2, lpcol1:lpcol2)= 1;%rand(sizefiltrow)*2*pi - pi; %debug

ImSwp(:,:,layer) = ...
ifft2(Amp(:,:,layer).*exp(sqrt(-1)*(ifftshift(rPhase(:,:,layer)) )));   
%combine Amp and Phase then perform inverse Fourier


ImSwp = real(ImSwp); %get rid of imaginary

imwrite(uint8(ImSwp), ...
    [newdirname imagename(1:4) '_' num2str(m) '.jpg'], ...
    'jpeg', 'quality', 100)

disp(['Mean intensity of original image is ' num2str(mean(mean(Img))) ])
disp(['Mean intensity of filtered image is ' num2str(mean(mean(fswap_Im))) ])
disp(['Mean intensity of phase-rand image is ' num2str(mean(mean(ImSwp))) ])


figure(1); clf;
hold on
subplot(2,3,1)

% imagesc(Img)
 imshow(uint8(Img))
% axis square 
title('Original image')

subplot(2,3,2)

% imagesc(fswap_Im)
 imshow(uint8(fswap_Im))
% axis square 
title(['Image filtered with ' num2str(lpcutoff)  'cpi ' upper(filtertype)])

subplot(2,3,3) 
%  imagesc(ImSwp)
 imshow(uint8(ImSwp))
% axis square
title(['Phase randomized for ' ...
   num2str(bpcutoff1) '-' num2str(bpcutoff2)  ' cpi '])
% set(gca,'xtick', [], 'ytick', []);
% figure(1)
% 
% figure(f2), 
subplot(2,3,4), 
 imagesc(fftshift(abs(log (1+ fft2( (Img),sizerow,sizecol) ) )))
 set(gca,'xtick', [], 'ytick', []);
title('Image FFT')
axis square

subplot(2,3,5), 
% imagesc(fftshift(abs(log (1+ fft2( (Img),sizerow,sizecol) ) )))
imshow(fftshift(Phase))
title('Original image phase')
axis square

subplot(2,3,6), 
imshow(rPhase)
title('Modified image phase')
figure(1)
%part in image (due to rounding error)
% imwrite(ImScrambled,'BearScrambled.jpg','jpg');

% imshow(ImScrambled)
% f2 = figure(2);
%   imshow(uint8(ImSwp))
% %   imshow(uint8(fswap_Im))
% % axis square
% % if strcmp('lp', filtertype)
% % title([upper(filtertype) '-filtered with ' ...
% %    num2str(lpcutoff)  ' cpi cutoff'])
% % elseif strcmp('hp', filtertype)
% %     title([upper(filtertype) '-filtered with ' ...
% %    num2str(hpcutoff)  ' cpi cutoff'])
% % end
% 
% % title(['Phase-randomized in the range ' ...
% %   num2str(bpcutoff1) '-' num2str(bpcutoff2)  ' cpi'])
% figure(f2)


% f3=figure(3);

% if strcmp('lp', filtertype)
% mesh(fftshift(L))
% title([num2str(lpcutoff) 'cpi ' filtername ' ' upper(filtertype)  ' filter'])
% elseif  strcmp('hp', filtertype)
%    mesh(fftshift(H))
%    title([num2str(hpcutoff) 'cpi ' filtername ' ' upper(filtertype)  ' filter'])
% end
% title([num2str(lpcutoff) 'cpi ' filtername ' ' upper(filtertype)  ' filter'])

% figure(f3)

end

end %m

  end %first i
  