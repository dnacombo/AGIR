%clear
% gimme your images
[listimages,pathname] = uigetfile({'*.jpg' '*.tif' '*.*'}','Select all image files to process','MultiSelect','on');
if isequal(listimages,0)
    return
elseif ischar(listimages)
    listimages = cellstr(listimages);
end
for i = 1:numel(listimages)
    listimages{i} = fullfile(pathname,listimages{i});
end

%%
for i_f = 1:length(listimages)% for each image
    try
        Im = double(imread([listimages{i_f}]));% read it
    catch
        disp(['Image ' [listimages{i_f}] ' unreadeable. Skipped']);
        continue
    end
    ImHSV = rgb2hsv(Im);
    % ImHSV(:,:,1) is Hue
    % ImHSV(:,:,2) is Saturation
    % ImHSV(:,:,3) is value
    
    % your contrast variables are here
    contrast = .8;%  0 --> gray
    %               .5 --> half contrast    
    %               1 --> normal
    %               >1 --> hyper contrast
    
    % change contrast
    Vs = ImHSV(:,:,3);
    ImHSV(:,:,3) = (ImHSV(:,:,3) - mean(Vs(:))) * contrast + mean(Vs(:));
    Im = uint8(hsv2rgb(ImHSV));
    
    imshow(Im);
    uiwait(gcf);
end


