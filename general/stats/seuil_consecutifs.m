function [data] = seuil_consecutifs(data,seuil,nb,dir)

% [data] = seuil_consecutifs(data,seuil,nb,dir)
%
% Retourne data(i) = abs(data(i)) si (abs(data(i)) > seuil) et 
% (nb datapoints autour de data(i) > seuil dans la 2è dimension)
% data(i) = 0 sinon
% dir donne le sens de la comparaison : '>' ou '<' (défaut = '>')
% dans l'explication ci-dessus, dir = '>'
%
% Max 31-10-2007

if length(size(data))>2
    error('data should have more than 2 dimensions');
end
if not(exist('dir','var'))
    dir = '>'; % dir
end

datawork = data;
% déjà on met tout ce qui dépasse pas le seuil (dans le bon sens) à NaN
switch dir
    case '>'
        datawork(find(abs(data) < seuil)) = NaN;
    case '<'
        datawork(find(abs(data) > seuil)) = NaN;
end
in = 0;
% ensuite, pour tout ce qui dépasse, on check que les nb echantillons
% suivants aussi.
% Sur chaque ligne, 
for i_sens = 1:size(data,1)
    i = 1;
    while i < size(data,2)-nb % pour toutes les colonnes
        while all(not(isnan(datawork(i_sens,i:i+nb)))) % tant que 
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
        datawork(i_sens,i) = NaN; % si pour cet échantillon, 
        % les nb suivants ne passent pas le seuil, on le met à NaN
        i = i+1;
    end
    datawork(i_sens,i:end) = NaN;% les derniers, on les met à NaN
end
data = datawork;
return



