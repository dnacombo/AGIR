function [Amp, Phase] = fourier_phase_amp(Im)

% [Amp, Phase] = myfourier(Im)
% retrieve Amplitude and phase of fourier transform of image.

Imsize = size(Im);
if numel(Imsize) == 2
    numLayers = 1;
elseif numel(Imsize) == 3
    numLayers = Imsize(3);
end

for i_layer = 1:numLayers
    %Fast-Fourier transform
    ImFourier(:,:,i_layer) = fft2(Im(:,:,i_layer));
    %amplitude spectrum
    Amp(:,:,i_layer) = abs(ImFourier(:,:,i_layer));
    %phase spectrum
    Phase(:,:,i_layer) = angle(ImFourier(:,:,i_layer));
end

