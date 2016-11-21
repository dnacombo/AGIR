function data = blend_transparency(data,mask,alpha,bcg)

% data = blend_transparency(data,mask,alpha,bcg)
% 
% data will be made transparent by a factor 1-alpha toward the color bcg
% only in points where mask == 0
% data is to be passed to surf before entering here and then to
% CDataMapping of the surface.
%
% Inputs :  data : NxM double data
%           mask : NxM double mask (zeros and ones)
%           alpha : color triplet for transparency of each RGB channel
%           bcg : color triplet of the background
%
% Output : data : NxMx3 color data to pass to surf...
%
% 

cmap = colormap;

alpha = repmat(alpha',[1,size(data)]);
alpha = permute(alpha,[2:length(size(alpha)),1]);
bcg = repmat(bcg',[1,size(data)]);
bcg = permute(bcg,[2:length(size(bcg)),1]);

maxdat = max(data(:));
mindat = min(data(:));
% data = (data - min(data(:)));
% data = data./max(data(:));
% data(find(data<0)) = 0;
% data(find(data>0)) = 1;
colorlimits = caxis;
% Datamap = linspace(colorlimits(1),colorlimits(2),size(cmap,1));

data = ceil((size(cmap,1)-1)/diff(colorlimits) .* data + (colorlimits(2) * 1 - size(cmap,1) * colorlimits(1)) / diff(colorlimits));

% data = (data + (mindat/colorlimits(1))) ./ (maxdat/colorlimits(2));
% 
% data = ceil(data.*(size(cmap,1)-1)) + 1;

data = reshape(cmap(data,:),[size(data),3]);

mask = repmat(mask,[1 1 3]);


%data = data + (1-mask) .* alpha;
id = find(1-mask);
data(id) = data(id) .* alpha(id) + bcg(id) .* (1 - alpha(id));







