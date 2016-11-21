function [TRs_out, varargout] = eject_outliers(TRs_in,TR_min,nsigma)

% TRs_out = eject_outliers(TRs_in,TR_min,nsigma)
% Nettoie TRs_in (double) en retirant les outliers.
% Retire les TRs de moins de TR_min (.3 par défaut)
% et les TRs en dehors de +/- nsigma (2 par défaut) ecarts type autour de la moyenne.
% 
% [TRs_out, idx] = eject_outliers(TRs_in,TR_min,nsigma)
% renvoie idx (double) l'index des TRs supprimés de TRs_in dans TRs_in.

if nargin < 2
    TR_min = .3;
end
if nargin < 3
    nsigma = 2;
end

TRs_work = TRs_in;

trouv1 = find(TRs_work < TR_min);% trouve les TRs trop courts
moy = nanmean(TRs_work);
ecarts = nanstd(TRs_work);
trouv2 = find(TRs_work > moy + nsigma*ecarts | TRs_work < moy - nsigma*ecarts);% trouve >< moy +- nsigma ecartypes

TRs_work([trouv1 trouv2]) = NaN; % efface les 2

TRs_out = TRs_work;

if nargout > 1
    varargout{1} = [trouv1 trouv2];
end
