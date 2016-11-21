% pop_eegfiltnew() - Filter data using Hamming windowed sinc FIR filter
%
% Usage:
%   >> [dataf, b] = eegfilt(data, param, value... );
%
% Inputs:
%   data       - EEG data (chans x time) or (chans x time x trials) 
% Param - value pairs
%   locutoff  - lower edge of the frequency pass band (Hz)
%               {[]/0 -> lowpass} 
%   hicutoff  - higher edge of the frequency pass band (Hz)
%               {[]/0 -> highpass}
%
% Optional inputs:
%   filtorder - filter order (filter length - 1). Mandatory even
%   revfilt   - [0|1] invert filter (from bandpass to notch filter)
%               {default 0 (bandpass)}
%   usefft    - ignored (backward compatibility only)
%   plotfreqz - [0|1] plot filter's frequency and phase response
%               {default 0} 
%   minphase  - scalar boolean minimum-phase converted causal filter
%               {default false}
%
% Outputs:
%   dataf     - filtered data
%   b         - filter coefficients
%
%
% Author: Andreas Widmann, University of Leipzig, 2012
%           adapted to not rely on EEGLAB by Max Chaumon 2016
%
% See also:
%   firfilt, firws, windows


% Copyright (C) 2008 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [dataf, b] = eegfilt(data, varargin)


def = [];
def.srate = [];
def.locutoff = [];
def.hicutoff = [];
def.filtorder = [];
def.revfilt = 0;
def.usefft = [];
def.plotfreqz = 0;
def.minphase = 0;

if all(cellfun(@isnumeric,varargin(1:3)))
    def.srate = varargin{1};
    def.locutoff = varargin{2};
    def.hicutoff = varargin{3};
    varargin(1:3) = [];
end
struct2ws(vararg2cfg(varargin,def,1));

% Constants
TRANSWIDTHRATIO = 0.25;
fNyquist = srate / 2;

% Check arguments
if locutoff == 0, locutoff = []; end
if hicutoff == 0, hicutoff = []; end
if isempty(hicutoff) % Convert highpass to inverted lowpass
    hicutoff = locutoff;
    locutoff = [];
    revfilt = ~revfilt;
end
edgeArray = sort([locutoff hicutoff]);

if isempty(edgeArray)
    error('Not enough input arguments.');
end
if any(edgeArray < 0 | edgeArray >= fNyquist)
    error('Cutoff frequency out of range');
end

if ~isempty(filtorder) && (filtorder < 2 || mod(filtorder, 2) ~= 0)
    error('Filter order must be a real, even, positive integer.')
end

% Max stop-band width
maxTBWArray = edgeArray; % Band-/highpass
if revfilt == 0 % Band-/lowpass
    maxTBWArray(end) = fNyquist - edgeArray(end);
elseif length(edgeArray) == 2 % Bandstop
    maxTBWArray = diff(edgeArray) / 2;
end
maxDf = min(maxTBWArray);

% Transition band width and filter order
if isempty(filtorder)

    % Default filter order heuristic
    if revfilt == 1 % Highpass and bandstop
        df = min([max([maxDf * TRANSWIDTHRATIO 2]) maxDf]);
    else % Lowpass and bandpass
        df = min([max([edgeArray(1) * TRANSWIDTHRATIO 2]) maxDf]);
    end

    filtorder = 3.3 / (df / srate); % Hamming window
    filtorder = ceil(filtorder / 2) * 2; % Filter order must be even.
    
else

    df = 3.3 / filtorder * srate; % Hamming window
    filtorderMin = ceil(3.3 ./ ((maxDf * 2) / srate) / 2) * 2;
    filtorderOpt = ceil(3.3 ./ (maxDf / srate) / 2) * 2;
    if filtorder < filtorderMin
        error('Filter order too low. Minimum required filter order is %d. For better results a minimum filter order of %d is recommended.', filtorderMin, filtorderOpt)
    elseif filtorder < filtorderOpt
        warning('firfilt:filterOrderLow', 'Transition band is wider than maximum stop-band width. For better results a minimum filter order of %d is recommended. Reported might deviate from effective -6dB cutoff frequency.', filtorderOpt)
    end

end

filterTypeArray = {'lowpass', 'bandpass'; 'highpass', 'bandstop (notch)'};
fprintf('eegfiltnew() - performing %d point %s filtering.\n', filtorder + 1, filterTypeArray{revfilt + 1, length(edgeArray)})
fprintf('eegfiltnew() - transition band width: %.4g Hz\n', df)
fprintf('eegfiltnew() - passband edge(s): %s Hz\n', mat2str(edgeArray))

% Passband edge to cutoff (transition band center; -6 dB)
dfArray = {df, [-df, df]; -df, [df, -df]};
cutoffArray = edgeArray + dfArray{revfilt + 1, length(edgeArray)} / 2;
fprintf('pop_eegfiltnew() - cutoff frequency(ies) (-6 dB): %s Hz\n', mat2str(cutoffArray))

% Window
winArray = windows('hamming', filtorder + 1);

% Filter coefficients
if revfilt == 1
    filterTypeArray = {'high', 'stop'};
    b = firws(filtorder, cutoffArray / fNyquist, filterTypeArray{length(cutoffArray)}, winArray);
else
    b = firws(filtorder, cutoffArray / fNyquist, winArray);
end

if minphase
    disp('eegfiltnew() - converting filter to minimum-phase (non-linear!)');
    b = minphaserceps(b);
end

% Plot frequency response
if plotfreqz
    freqz(b, 1, 8192, srate);
end

% Filter
if minphase
    disp('pop_eegfiltnew() - filtering the data (causal)');
    dataf = firfiltsplit(data, b, 1);
else
    disp('pop_eegfiltnew() - filtering the data (zero-phase)');
    dataf = firfilt(data, b);
end


% firfiltsplit() - Split data at discontinuities and forward to dc padded
%                  filter function
%
% Usage:
%   >> dataf = firfiltsplit(data, b);
%
% Inputs:
%   data           - EEG data
%   b             - vector of filter coefficients
%   causal        - scalar boolean perform causal filtering {default 0}
%
% Outputs:
%   dataf           - filtered EEG data
%
% Note:
%   This function is (in combination with firfiltdcpadded) just a
%   non-memory optimized version of the firfilt function allowing causal
%   filtering. Will possibly replace firfilt in the future.
%
% Author: Andreas Widmann, University of Leipzig, 2013
%
% See also:
%   firfiltdcpadded, findboundaries

% Copyright (C) 2013 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function data = firfiltsplit(data, b, causal)

if nargin < 3 || isempty(causal)
    causal = 0;
end
if nargin < 2
    error('Not enough input arguments.');
end

s = size(data);
% Find data discontinuities and reshape epoched data
if numel(s) == 3 % Epoched data
    data = data(:,:);
    dcArray = 1 : s(2) : s(2) * (s(3) + 1);
else % Continuous data
    dcArray = [1 s(2) + 1];
end

% Loop over continuous segments
for iDc = 1:(length(dcArray) - 1)

    % Filter segment
    data(:, dcArray(iDc):dcArray(iDc + 1) - 1) = firfiltdcpadded(b, data(:, dcArray(iDc):dcArray(iDc + 1) - 1)', causal)';

end

% Reshape epoched data
data = reshape(data, s);

% firfiltdcpadded() - Pad data with DC constant and filter
%
% Usage:
%   >> data = firfiltdcpadded(data, b, causal);
%
% Inputs:
%   data      - raw data
%   b         - vector of filter coefficients
%   causal    - boolean perform causal filtering {default 0}
%
% Outputs:
%   data      - smoothed data
%
% Note:
%   firfiltdcpadded always operates (pads, filters) along first dimension.
%   Not memory optimized.
%
% Author: Andreas Widmann, University of Leipzig, 2013
%
% See also:
%   firfiltsplit

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2013 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [ data ] = firfiltdcpadded(b, data, causal)

% Defaults
if nargin < 3 || isempty(causal)
    causal = 0;
end

% Check arguments
if nargin < 2
    error('Not enough input arguments.');
end

% Filter's group delay
if mod(length(b), 2) ~= 1
    error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;
b = double(b); % Filter with double precision

% Pad data with DC constant
if causal
    startPad = repmat(data(1, :), [2 * groupDelay 1]);
    endPad = [];
else
    startPad = repmat(data(1, :), [groupDelay 1]);
    endPad = repmat(data(end, :), [groupDelay 1]);
end

% Filter data
data = filter(b, 1, double([startPad; data; endPad])); % Pad and filter with double precision

% Remove padded data
data = data(2 * groupDelay + 1:end, :);
 

% firfilt() - Pad data with DC constant, filter data with FIR filter,
%             and shift data by the filter's group delay
%
% Usage:
%   >> dataf = firfilt(data, b, nFrames);
%
% Inputs:
%   data          - EEG data
%   b             - vector of filter coefficients
%
% Optional inputs:
%   nFrames       - number of frames to filter per block {default 1000}
%
% Outputs:
%   dataf         - filtered EEG data
%
% Note:
%   Higher values for nFrames increase speed and working memory
%   requirements.
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   filter, findboundaries

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function data = firfilt(data, b, nFrames)

if nargin < 2
    error('Not enough input arguments.');
end
if nargin < 3 || isempty(nFrames)
    nFrames = 1000;
end

% Filter's group delay
if mod(length(b), 2) ~= 1
    error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;

s = size(data);
% Find data discontinuities and reshape epoched data
if numel(s) == 3 % Epoched data
    data = data(:,:);
    dcArray = 1 : s(2) : s(2) * (s(3) + 1);
else % Continuous data
    dcArray = [1 s(2) + 1];
end

% Initialize progress indicator
nSteps = 20;
step = 0;
fprintf(1, 'firfilt(): |');
strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
tic

for iDc = 1:(length(dcArray) - 1)

        % Pad beginning of data with DC constant and get initial conditions
        ziDataDur = min(groupDelay, dcArray(iDc + 1) - dcArray(iDc));
        [temp, zi] = filter(b, 1, double([data(:, ones(1, groupDelay) * dcArray(iDc)) ...
                                  data(:, dcArray(iDc):(dcArray(iDc) + ziDataDur - 1))]), [], 2);

        blockArray = [(dcArray(iDc) + groupDelay):nFrames:(dcArray(iDc + 1) - 1) dcArray(iDc + 1)];
        for iBlock = 1:(length(blockArray) - 1)

            % Filter the data
            [data(:, (blockArray(iBlock) - groupDelay):(blockArray(iBlock + 1) - groupDelay - 1)), zi] = ...
                filter(b, 1, double(data(:, blockArray(iBlock):(blockArray(iBlock + 1) - 1))), zi, 2);

            % Update progress indicator
            [step, strLength] = mywaitbar((blockArray(iBlock + 1) - groupDelay - 1), size(data, 2), step, nSteps, strLength);
        end

        % Pad end of data with DC constant
        temp = filter(b, 1, double(data(:, ones(1, groupDelay) * (dcArray(iDc + 1) - 1))), zi, 2);
        data(:, (dcArray(iDc + 1) - ziDataDur):(dcArray(iDc + 1) - 1)) = ...
            temp(:, (end - ziDataDur + 1):end);

        % Update progress indicator
        [step, strLength] = mywaitbar((dcArray(iDc + 1) - 1), size(data, 2), step, nSteps, strLength);

end

% Reshape epoched data
data = reshape(data, s);

% Deinitialize progress indicator
fprintf(1, '\n')


function [step, strLength] = mywaitbar(compl, total, step, nSteps, strLength)

progStrArray = '/-\|';
tmp = floor(compl / total * nSteps);
if tmp > step
    fprintf(1, [repmat('\b', 1, strLength) '%s'], repmat('=', 1, tmp - step))
    step = tmp;
    ete = ceil(toc / step * (nSteps - step));
    strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '%s %3d%%, ETE %02d:%02d'], progStrArray(mod(step - 1, 4) + 1), floor(step * 100 / nSteps), floor(ete / 60), mod(ete, 60));
end


