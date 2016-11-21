function [ results ] = mylimo_contrast( Y, Betas, LIMO, TF, analysistype, targetcontrast)
%mylimo_contrast
%   what it says.

if not(exist('targetcontrast','var'))
    targetcontrast = size(LIMO.contrast,2);
end

switch analysistype
    case 1
        
        X = LIMO.design.X;
        
        C = LIMO.contrast{targetcontrast}.C;
        % compute Projection onto the error
        R = eye(size(Y,1)) - (X*pinv(X));
        
        dfe     = LIMO.model.model_df(2);
        switch TF
            case 'T'
                
                var   = ((R*Y)'*(R*Y)) / dfe;
                
                con(:,1) = C*Betas ;  % contrast
                con(:,2) = (C*Betas) ./ sqrt(diag(var)'.*(C*pinv(X'*X)*C')); % T value
                con(:,3) = 2*(1-tcdf(squeeze(con(:,2)), dfe)); % p value
                
                results = con;
                
            case 'F'
                ess(:,1:size(C,1)) = (C*Betas)' ; % contrast
                
                E = (Y'*R*Y);
                c = zeros(length(C));
                c(1:size(C,1),1:size(C,2)) = C;
%                 for n=1:length(C)
%                     c(n,n) = C(n,n);
%                 end
                
                try
                    C0 = eye(rank(X)+1) - C*pinv(C);
                catch ME
                    C0 = eye(rank(X)) - c*pinv(c);
                end
                X0 = X*C0;
                R0 = eye(size(Y,1)) - (X0*pinv(X0));
                M = R0 - R;
                H = (Betas'*X'*M*X*Betas);
                if rank(c) == 1
                    df = 1;
                else
                    df = rank(c) - 1;
                end
                ess(:,end-1)    = (diag(H)/df)./(diag(E)/dfe);  % F value
                ess(:,end) = 1 - fcdf(ess(:,end-1), rank(c)-1, dfe);   % p value
                
                results = ess;
                
        end
        
    otherwise
        error('not implemented')
end


