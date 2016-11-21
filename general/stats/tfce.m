function tfce_score = tfce(data,varargin)

% Implementation of the Threshold-free cluster enhancement method
% developped for fMRI by Smith & Nichols, NeuroImage 44(2009), 83-98
%
% INPUT tfce_score = tfce(data)
%       tfce_score = tfce(data,channeighbstructmat)
%       tfce_score = tfce(data,channeighbstructmat,dh,E,H)
%
%       data is NDimensional map of values
%       if data has more than one element along the first dimension
%       i.e. size(data,1) > 1, this dimension can take connectivity
%       described in channeighbstructmat (if this is empty, this feature is
%       ignored.)
%       dh, E, H and are the parameters of the tfce algorithm defaults are 0.1, 0.5, 2
%       tfce = sum(extent(h)^E*height^H*dh)
%
% OUPUT tfce_score is a map of scores
%
%

% Maximilien Chaumon
% Based on Cyril Pernet 18-10-2011


max_precision = 200; % define how many thresholds between min t/F map and
% max t/F map --> needed as sometime under H0 some values can be
% arbritrarily high due to low variance during resampling

%% check input

if nargin > 1
    channeighbstructmat = varargin{1};
else
    channeighbstructmat = [];
end
if nargin > 2
    dh = varargin{2};
else
    dh = .1;
end
if nargin > 3
    E = varargin{3};
else
    E = .5;
end
if nargin > 4
    H = varargin{4};
else
    H = 2;
end
if ~isempty(channeighbstructmat) && ...
        (size(channeighbstructmat,1) ~= size(channeighbstructmat,2)...
        || size(channeighbstructmat,1) ~= size(data,1))
    error('Wrong size of channeighbstructmat')
end

% check negative values if so do negate and add scores
if min(data(:)) < 0
    neg_data = -(data < 0) .* data;
    neg_tfce = tfce(neg_data,channeighbstructmat,dh,E,H);
    pos_data = (data > 0) .* data;
    pos_tfce = tfce(pos_data,channeighbstructmat,dh,E,H);
    tfce_score = neg_tfce + pos_tfce;
    return
end

s = size(data);


% start tcfe

% define increment size forced by dh
data_range  = range(data(:));
precision = round(data_range / dh);
% if precision > max_precision % arbitrary decision to limit precision to 200th of the data range - needed as sometime under H0 one value can be very wrong
%     precision = max_precision;
% end
increment = data_range / precision;

% select a height, obtain cluster map, obtain extent map (=cluster
% map but with extent of cluster rather than number of the cluster)
% then tfce score for that height
steps = min(data(:)):increment:max(data(:));
nsteps = numel(steps);
% inserting steps as first dimension.
tfce_score = NaN([nsteps s]);
% f = waitbar(0,'Thresholding levels','name','TFCE');
for isteps = 1:nsteps
%     waitbar(isteps/nsteps);
    if s(1) == 1
        [C,nb] = bwlabeln((data > steps(isteps)));
    else
        [C nb] = findcluster(data > steps(isteps), channeighbstructmat);
    end
    extent_map = zeros(size(C)); % same as cluster map but contains extent value instead
    for i = 1:nb
        idx = C(:) == i;
        extent_map(idx) = extent_map(idx) + sum(idx);
    end
    tfce_score(isteps,:) = (extent_map(:).^E) .* steps(isteps) .^ H .* dh;
end
% compute final score
% repermute to suppress first dimension (was nsteps)
tfce_score = permute(nansum(tfce_score,1),[2:numel(size(tfce_score)) 1]);
% close(f)

function [cluster, num] = findcluster(onoff, spatdimneighbstructmat, varargin)

% FINDCLUSTER returns all connected clusters in a 3 dimensional matrix
% with a connectivity of 6.
%
% Use as
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, minnbchan)
% or as
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, spatdimneighbselmat, minnbchan)
% where
%   onoff                   is a 3D boolean matrix with size N1xN2xN3,
%                           N1=number of channels
%                           N2 & N3 = Time Frequency (any order)
%   spatdimneighbstructmat  defines the neighbouring channels/combinations, see below
%   minnbchan               the minimum number of neighbouring channels/combinations
%   spatdimneighbselmat     is a special neighbourhood matrix that is used for selecting
%                           channels/combinations on the basis of the minnbchan criterium
%
% The neighbourhood structure for the first dimension is specified using
% spatdimneighbstructmat, which is a 2D (N1xN1) matrix. Each row and each column corresponds
% to a channel (combination) along the first dimension and along that row/column, elements
% with "1" define the neighbouring channel(s) (combinations). The first dimension of
% onoff should correspond to the channel(s) (combinations).
% The lower triangle of spatdimneighbstructmat, including the diagonal, is
% assumed to be zero.
%
% See also BWSELECT, BWLABELN (image processing toolbox)
% and SPM_CLUSTERS (spm2 toolbox).

% Copyright (C) 2004, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: findcluster.m 952 2010-04-21 18:29:51Z roboos $
% renamed for integration in LIMO toolbox: GAR, University of Glasgow, June
% 2010

spatdimlength = size(onoff, 1);
nfreq = size(onoff, 2);
ntime = size(onoff, 3);

if length(size(spatdimneighbstructmat))~=2 || ~all(size(spatdimneighbstructmat)==spatdimlength)
    error('invalid dimension of spatdimneighbstructmat');
end

minnbchan=0;
if length(varargin)==1
    minnbchan=varargin{1};
end;
if length(varargin)==2
    spatdimneighbselmat=varargin{1};
    minnbchan=varargin{2};
end;

if minnbchan>0
    % For every (time,frequency)-element, it is calculated how many significant
    % neighbours this channel has. If a significant channel has less than minnbchan
    % significant neighbours, then this channel is removed from onoff.
    
    if length(varargin)==1
        selectmat = single(spatdimneighbstructmat | spatdimneighbstructmat');
    end;
    if length(varargin)==2
        selectmat = single(spatdimneighbselmat | spatdimneighbselmat');
    end;
    nremoved=1;
    while nremoved>0
        nsigneighb=reshape(selectmat*reshape(single(onoff),[spatdimlength (nfreq*ntime)]),[spatdimlength nfreq ntime]);
        remove=(onoff.*nsigneighb)<minnbchan;
        nremoved=length(find(remove.*onoff));
        onoff(remove)=0;
    end;
end;

% for each channel (combination), find the connected time-frequency clusters
labelmat = zeros(size(onoff));
total = 0;
for spatdimlev=1:spatdimlength
    [labelmat(spatdimlev, :, :), num] = bwlabeln(reshape(onoff(spatdimlev, :, :), nfreq, ntime), 4);
    labelmat(spatdimlev, :, :) = labelmat(spatdimlev, :, :) + (labelmat(spatdimlev, :, :)~=0)*total;
    total = total + num;
end

% combine the time and frequency dimension for simplicity
labelmat = reshape(labelmat, spatdimlength, nfreq*ntime);

% combine clusters that are connected in neighbouring channel(s)
% (combinations).
replaceby=1:total;
for spatdimlev=1:spatdimlength
    neighbours=find(spatdimneighbstructmat(spatdimlev,:));
    for nbindx=neighbours
        indx = find((labelmat(spatdimlev,:)~=0) & (labelmat(nbindx,:)~=0));
        for i=1:length(indx)
            a = labelmat(spatdimlev, indx(i));
            b = labelmat(nbindx, indx(i));
            if replaceby(a)==replaceby(b)
                % do nothing
                continue;
            elseif replaceby(a)<replaceby(b)
                % replace all entries with content replaceby(b) by replaceby(a).
                replaceby(replaceby==replaceby(b)) = replaceby(a);
            elseif replaceby(b)<replaceby(a)
                % replace all entries with content replaceby(a) by replaceby(b).
                replaceby(replaceby==replaceby(a)) = replaceby(b);
            end
        end
    end
end

% renumber all the clusters
num = 0;
cluster = zeros(size(labelmat));
for uniquelabel=unique(replaceby(:))'
    num = num+1;
    cluster(ismember(labelmat(:),find(replaceby==uniquelabel))) = num;
end

% reshape the output to the original format of the data
cluster = reshape(cluster, spatdimlength, nfreq, ntime);
