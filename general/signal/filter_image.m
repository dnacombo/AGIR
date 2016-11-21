function Imf = filter_image(Im, type, cutoff, rescale)
% Imf = filter_image(Im, type, cutoff)
% 'low' or 'high' pass filter image at cutoff cpi using a gaussian filter.

if not(exist('rescale', 'var'))
    rescale = 1;
end

if length(size(Im)) > 2
    error('Image data should be grayscale. (only 2 dimensions)')
end

picsize=size(Im);


Im = double(Im);

switch type
        
    case 'low'

        lp_cutoff = round(cutoff);
        L = lpfilter('gaussian',picsize(1),picsize(2),lp_cutoff,1);
        %           L = lpfilter('btw',picsize(1),picsize(2),lp_cutoff,1);
        Imf = dftfilt(Im,L);
%         mean(mean(Imf))
%         max(max(Imf))
%         min(min(Imf))
        %set intensity to mean
%         Imf =Imf-Immean;
    case  'high'
        hp_cutoff = round(cutoff);
        H = hpfilter('gaussian',picsize(1),picsize(2),hp_cutoff,1);
        %             H = hpfilter('btw',picsize(1),picsize(2),hp_cutoff,1);
%         Imf = dftfilt(Im,H);
%         Imf = Imf+Immean;
end

if rescale == 2000
    % rescale so that values of Imf have same max and min as Im
    Immean = mean(mean(Im));
    Immax = max(max(Im));
    Immin = min(min(Im));

    
    Imfmax = max(max(Imf));
    Imfmin = min(min(Imf));

    Imf = (Imfmax - Imfmin)/(Immax - Immin) * Imf + Imfmin;
    Imfmean =  mean(mean(Imf));
    
%     Imfmean == Immean
end