function Ap = pHF_to_Aprime(pHF)

% Ap = pHF_to_Aprime(pHF)
% converts proportion Hit and FA to A-prime

H = pHF(:,1);
F = pHF(:,2);

% formula from equation (3) in 
% Stanislaw H, Todorov N (1999) Calculation of signal detection theory measures. 
% Behavior Research Methods, Instruments, & Computers 31:137–149.
Ap = .5 + (sign(H-F) .* ((H-F).^2 + abs(H-F))./(4.*max(H,F) - 4.*H.*F));

