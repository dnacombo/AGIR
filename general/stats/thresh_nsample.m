function [mask] = thresh_nsample(data,thresh,nb,comp)

% mask = thresh_nsample(data,thresh,nb,comp)
%
% Returns a mask the same size as data which is true where
% eval([data comp seuil]) is true for at least nb datapoints along
% dimension 2. 
% comp is the comparaison : '>' '<' '==' etc. (default = '<')
%
% input: data = N x M matrix of data to compare
%       thresh = threshold of comparison should be either a scalar or the
%               same size as data
%       nb = nb of samples for which comparison must be true along
%               dimension 2
%       comp = the comparison e.g. '>' '<' '=='...
%
% Max 31-10-2007
% v2: 2010-12

if ndims(data)>2
    error('data should have no more than 2 dimensions');
end
if not(exist('comp','var'))
    comp = '<';
end

mask = eval(['data' comp 'thresh']);

in = 0;
% ensuite, pour tout ce qui dépasse, on check que les nb echantillons
% suivants aussi.
% Sur chaque ligne, 
for i_sens = 1:size(data,1)
    i = 1;
    while i < size(data,2)-nb % pour toutes les colonnes
        while all(mask(i_sens,i:i+nb)) % tant que 
            % tous les nb échantillons suivants dépassent le seuil
            i = i+1;% on regarde le suivant
            in = 1; % quoi qu'il arrive on laissera nb échantillons à leur valeur initiale
            if i == size(data,2)-nb
                break
            end
        end
        if in
            i = i+nb; % on laisse nb échantillons à leur valeur initiale
            in = 0;
        end
        mask(i_sens,i) = 0; % si pour cet échantillon, 
        % les nb suivants ne passent pas le seuil, on le met à 0
        i = i+1;
    end
    mask(i_sens,i:end) = 0;% les derniers, on les met à 0
end
return



