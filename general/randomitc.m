function itcs = randomitc(ntri,nboot)

%itcs = randomitc(ntri,nboot)
% retrieve expeted itc distribution when no phase coherence;

if nargin == 0
    ntri = 100;
    nboot = 100;
end

phases = 2*pi*rand(nboot,ntri,'single');

itcs = single(abs(mean(exp(1i*phases),2)));



