function tmpeeglab = myresample(data, pnts, new_pnts);

if length(data) < 2
    tmpeeglab = data;
    return;
end;

%if size(data,2) == 1, data = data'; end;
% padding to avoid artifacts at the beginning and at the end
% Andreas Widmann May 5, 2011

%The pop_resample command introduces substantial artifacts at beginning and end
%of data when raw data show DC offset (e.g. as in DC recorded continuous files)
%when MATLAB Signal Processing Toolbox is present (and MATLAB resample.m command
%is used).
%Even if this artifact is short, it is a filtered DC offset and will be carried
%into data, e.g. by later highpass filtering to a substantial amount (easily up
%to several seconds).
%The problem can be solved by padding the data at beginning and end by a DC
%constant before resampling.

[p, q] = rat(pnts / new_pnts, 1e-12); % Same precision as in resample
N = 10; % Resample default
nPad = ceil((max(p, q) * N) / q) * q; % # datapoints to pad, round to integer multiple of q for unpadding
tmpeeglab = resample([data(ones(1, nPad), :); data; data(end * ones(1, nPad), :)], pnts, new_pnts);
nPad = nPad * p / q; % # datapoints to unpad
tmpeeglab = tmpeeglab(nPad + 1:end - nPad, :); % Remove padded data

