function [amps, freqs] = my_fft(data, timedim, srate, npoints);
%function [amps, freqs] = my_fft(data, timedim, srate, npoints)
%
% Function for computing FFT on EEG data. The function return only the
% positive half of the FFT spectrum - this is the part that you are
% normally interested in.
%
% - data: EEG data with arbitrary dimensions. Can be channels x times x
% trials or just a single vector of time points or whatever.
% - timedim: which dimension corresponds to time.
% - srate: sampling rate.
% npoints: length of the FFT
%
% Example:
% [amps, freqs] = my_fft(EEG.data, 2, EEG.srate, 1024)
%
% Written by Niko Busch - Charité Berlin (niko.busch@gmail.com)
%
% 2011-05-23


% Create the frequency axis.
if mod(npoints,2)==0
    k=-npoints/2:npoints/2-1; % N even
else
    k=-(npoints-1)/2:(npoints-1)/2; % N odd
end
T=npoints/srate;
freqs=k/T; 


% Compute and normalize FFT.
X = fft(data, npoints, timedim)/size(data,timedim); 


% Extract only the positive half of the spectrum
cutOff = ceil(npoints/2);

Xstring = 'X = X(';
for d = 1:ndims(X)
    if d==timedim
        Xstring = [Xstring '1:cutOff'];
    else
        Xstring = [Xstring ':'];
    end
    
    if d<ndims(X)
        Xstring = [Xstring ','];
    end
end
Xstring = [Xstring ');'];
eval(Xstring);

freqs = freqs(cutOff+1:end);


% Compute the power.
amps = abs(X);


