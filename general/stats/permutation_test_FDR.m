function [thersh_FWER,T_orig_ori,pthres,thresh_FDR_BH,thresh_FDR_BY] = permutation_test_FDR(data1,data2,M,label)

% [thersh_FWER,T_orig_ori,pthres,thresh_FDR_BH,thresh_FDR_BY] =
% permutation_test_FDR(data1,data2,M,label);
%
% permutation test --> data1 # data2
% permutation along dimension #1
%
% one input parameter data1 where are the differences between two conditions:
% H0 :  test for matched-pair = 0
% H1 :  data # 0 
%
% two input parameters data1, data2 
% case 1
%   H0: data1 =  data2  two independant samples and no differencial effect.
%       data1 and data2 are from the same population
%   H1: data1 # data2 and have diffent patterns
%
% three input parameters data1, data2, label 
% case 2  with  label 'matched-pair'   
%   H0:  data1-data2 two related samples and no differencial effect # 0
%       on the two observed behaviors
%   H1: data1 # data2 and have diffent patterns
%
% statistical measure T_orig = abs(Student value)
%
%input
% data1 : data condition A
%    if only one array for input parameter -> repeated measurements 
% data2 : data condition B
% M  number of resamples
% label : 'matched-pair'
%
%
%
%                       =BY
%                       Benjamini & Yekutieli (2001),"Resampling-based false discovery rate
%                       controlling multiple test proceedures for correleated test
%                       statistics".  J of Statistical Planning and Inference, 82:171-196.%
% 
% if label is omitting or missing label = independant 
%
% output 
% T_orig   statistical mapping for the difference between A and B
% thersh_FWER  FWER threshold computed for an alpha 
%
%  threh_FDR_BY threshold computed for False Discovery Rate  proportion of errors
%  FDR_option    = BH 
 %                        bt default Benjamini & Hochberg (1995), "Controlling the False Discovery Rate: A
%                        Practical and Powerful Approach to Multiple Testing". J Royal Stat Soc,
%                       Ser B.  57:289-300.
%                       The Benjamini & Hochberg (1995) False Discovery Rate (FDR) procedure
%                       finds a threshold u such that the expected FDR is at most q.  spm_P_FDR
%
%  threh_FDR_BY threshold computed for False Discovery Rate  proportion of
%  errors
%                       =BY
%                       Benjamini & Yekutieli (2001),"Resampling-based false discovery rate
%                       controlling multiple test proceedures for correleated test
%                       statistics".  J of Statistical Planning and Inference, 82:171-196.%
%
%
%
% Background
%
% For a given threshold on a statistic image, the False Discovery Rate
% is the proportion of suprathreshold voxels which are false positives.
% Recall that the thresholding of each voxel consists of a hypothesis
% test, where the null hypothesis is rejected if the statistic is larger
% than threshold.  In this terminology, the FDR is the proportion of
% rejected tests where the null hypothesis is actually true.
%
% A FDR procedure produces a threshold that controls the expected FDR
% at or below q.  The FDR adjusted p-value for a variable (pixel) is the smallest q
% such that the pixel would be suprathreshold.
%
% In comparison, a traditional multiple comparisons proceedure
% (e.g. Bonferroni or random field methods) controls Familywise Error
% rate (FWER) at or below alpha.  FWER is the chance of one or more
% false positives *anywhere* (not just among suprathreshold pixels).  A
% FWER adjusted p-value for a pixel is the smallest alpha such that the
% voxel would be suprathreshold.
%
% 
% If there is truely no signal in the image anywhere, then a FDR
% procedure controls FWER, just as Bonferroni and random field methods
% do. (Precisely, controlling E(FDR) yields weak control of FWE).  If
% there *is* some signal in the image, a FDR method will be more powerful
% than a traditional method.
%
%
% For careful definition of FDR-adjusted p-values (and distinction between
% corrected and adjusted p-values) see Yekutieli & Benjamini (1999).
%
%
%
% References
%
% Benjamini & Hochberg (1995), "Controlling the False Discovery Rate: A
% Practical and Powerful Approach to Multiple Testing". J Royal Stat Soc,
% Ser B.  57:289-300.
%
% Benjamini & Yekutieli (2001), "The Control of the false discovery rate
% in multiple testing under dependency". To appear, Annals of Statistics.
% Available at http://www.math.tau.ac.il/~benja 
%
% Yekutieli & Benjamini (1999). "Resampling-based false discovery rate
% controlling multiple test proceedures for correleated test
% statistics".  J of Statistical Planning and Inference, 82:171-196.
%____________________________________________________________________
% /---Script Authors--------------------------------------\
% |                                                      |
% | *** Sylvain Baillet Ph.D.; Karim Jerbi               |
% | Cognitive Neuroscience & Brain Imaging Laboratory    |
% | CNRS UPR640 - LENA                                   | 
% | Hopital de la Salpetriere, Paris, France             |
% | sylvain.baillet@chups.jussieu.fr    
% | modified  by Jacques Martinerie   14/03/2005  25/12/2006
%    CNRS LENA |
% |                                                      |
% \------------------------------------------------------/
% %size(data1,1); %ntrials condition A
%size(data2,1); %ntrials condition B
%size(data1,2); % Number of sources/channel
%size(data1,3); %% Number dimensions to include in the test (time or frequency bins...)

% Date of creation: December 2004
% Date of revision version: Dec 2006

%initialize parameters
pthres = [0.05, 0.01, 0.005, 0.001];  % p-value thresholds selected
Nbmax=12;  %  maximum number of subjects for Fisher statistic

%M = 500; %nResample
nb_par=size(size(data1),2);
if size(data1,nb_par) == 1; nb_par=nb_par-1,  end;
JA = size(data1,1); %ntrials condition A
K = size(data1,2);  % Number of sources/channel
if exist('M') == 0;  M=500; end
if exist('label') == 0; label='unknown     '; end
if exist('FDR_option') == 0;  FDR_option='BH'; end

if (label ~= 'matched-pair');
    if exist('data2') == 1;
        JB = size(data2,1); %ntrials condition B
        if (JA+JB) <= Nbmax ;M=factorial(JA+JB)/(factorial(JA)*factorial(JB))-1;  end
        
    else
        label='matched-pair';
        data2=zeros(size(data1));
        if (size(data1,1)) < Nbmax ; M=2^(size(data1,1))-1;    ;  end % p-value thresholds selected 
        
    end
else
    if (JA) <= Nbmax ; M=2^(size(data1,1))-1; end
end
%if nb_par >= 3; F = size(data1,3);end %% Number dimensions to include in the test (time or frequency bins...)
%if nb_par == 4; TT=size(data1,4);end 
%if nb_par > 4; 'size Matrix  too large',return ;end 
if (floor(1/pthres(1))) >  M  ; pthres(1)=1/M; end



Siz=ndims(data1);
N_data1(1:Siz)=size(data1);
Narray=1; for k=2:Siz; Narray=Narray*N_data1(k); end; 
data1=reshape(data1,N_data1(1),Narray);
Siz=ndims(data2);
N_data2(1:Siz)=size(data2);
Narray=1; for k=2:Siz; Narray=Narray*N_data2(k); end; 
data2=reshape(data2,N_data2(1),Narray);
K=size(data2,2);   %  number of variables
gtot(1:JA)=0;
if strcmp(label,'matched-pair') ~= 1;      
    
    PB_orig(1:K) = mean(data2);
    PA_orig(1:K) = mean(data1);
    
    std_PA_orig(1:K)= std(data1);
    std_PB_orig(1:K)= std(data2);
    
    
    STD_PA_PB_orig=sqrt((std_PA_orig.*std_PA_orig/JA)+(std_PB_orig.*std_PB_orig/JB));
    % if abs( PA_orig - PB_orig ) > 1.e-06;     %  petite valeurs on néglige
    T_orig=( PA_orig - PB_orig ) ./ (STD_PA_PB_orig);
    % else
    %     T_orig=( PA_orig - PB_orig );
    %  end
    T_orig_ori=T_orig;
    T_orig=abs(T_orig);
else
    %matched-pairs
    data=data1-data2;
    
    
    PA_orig(1:K) = mean(data);
    std_PA_orig(1:K)= std(data);  
    
    STD_PA_PB_orig=sqrt(std_PA_orig.*std_PA_orig/JA);
    % if abs( PA_orig) > 1.e-06;     %  petite valeurs on néglige
    T_orig=(PA_orig)  ./ (STD_PA_PB_orig);
    % else
    % T_orig=(PA_orig);
    % end
    T_orig_ori=T_orig;
    T_orig=abs(T_orig);
end
% This is the statistics you want to test whether it departs from 0
Nb_H0=zeros(size(T_orig));
if M==0 ; thersh_FWER(1:4)=0;return; end

%Create Permutation Samples---------------------------------------------------

%save('data/with_pvalue/PermMatrix','PermMatrixA','PermMatrixB');
disp('Create permutation matrices -> Done')

%Create Perm Statistics---------------------------------------------------------
disp('Creating permutation statistics...')
hwait= waitbar(0,'Creating permutation statistics...') ;
Ps= T_orig;
data=[data1; data2];
for m = 1:M %for all permutations
    clear  PB_perm  PA_perm  std_PA_perm std_PB_perm
    if(~rem(m,10))
        waitbar(m/M,hwait)
        pack
    end
    j=randperm(size(data,1));
    if strcmp(label,'matched-pair') ~= 1; 
        if (JA+JB) <= Nbmax ;   %  number of subjects for random draft
            if m == 1;           
                n1=JA;
                n2=JB;
                nk=n1+n2;
                ij=1;
                j=2^nk-1;
                clear c1 c3
                for k=1:j
                    c= dec2bin(k);
                    a=size(c,2);
                    
                    c3(1:nk)=0;
                    ic=0;
                    for m1=a:-1:1
                        
                        c2=str2num(c(m1));
                        ic=ic+1;
                        
                        
                        c3(nk-ic+1)=c2;
                        if c2 > 0;
                            
                            c11=c3;
                            
                        end
                        
                        
                    end 
                    if sum(c11) == JA; c1(ij,:)=c11; ij=ij+1;end
                    
                end
                
                %  l1=find(c1 == 1);  l2=find(c1 == 0);
            end
            l1=find(c1(m,:) == 1);  l2=find(c1(m,:) == 0);
            
            
            data1_perm=data(l1,:);
            data2_perm=data(l2,:); 
            
            PB_perm(1:K) = mean(data2_perm);
            PA_perm(1:K) = mean(data1_perm);
            
            std_PA_perm(1:K)= std(data1_perm);
            std_PB_perm(1:K)= std(data2_perm);       
            
        else
            
            
            data1_perm=data(j(1:JA),:);
            data2_perm=data(j(JA+1:end),:); 
            
            PB_perm(1:K) = mean(data2_perm);
            PA_perm(1:K) = mean(data1_perm);
            
            std_PA_perm(1:K)= std(data1_perm);
            std_PB_perm(1:K)= std(data2_perm);        
            
        end
        STD_PA_PB_perm=sqrt((std_PA_perm.*std_PA_perm./JA)+(std_PB_perm.*std_PB_perm/JB));
        
        %   if abs( PA_perm - PB_perm ) > 1.e-06;     %  petite valeurs on néglige
       % T{m} =abs( PA_perm - PB_perm ) ./ (STD_PA_PB_perm);
          T{1} =abs( PA_perm - PB_perm ) ./ (STD_PA_PB_perm);
        %  else
        % T{m} =abs( PA_perm - PB_perm ) ;
        %end
    else
        % MATCHED-PAIR
        
        if exist('data2') ~ 0;
            if size(data1,1) > Nbmax ;   %  number of subjects for random draft
                
                clear g
                j=randperm(size(data1,1));
                %j1=randperm(size(data1,1));
                j1=floor(rand*(JA-1)+1);
                lm=j(1:j1(1));
                %if rem(m,2) ;lm=j(1:j1(1)); else; lm=j(end:-1:end-j1(end)+1); end
                %lm=find (j(1:JA) > JA/2); 
                g(1:JA)=1;
                g(lm)=-1;
                gtot=gtot+g;
                g=g';
                j=randperm(size(data1,1));
                %data1=data1(j);
                %data2=data2(j);  %  pour compenser les effets du tirage sous H0 
                data1_perm=data1-data2;   
                %gtot=gtot+data1_perm';
                
                for k=1:K
                    data1_perm(:,k)=data1_perm(:,k).*g;   
                end
                PA_perm(1:K) = mean(data1_perm);
                std_PA_perm(1:K)= std(data1_perm);     
                
                
            else
                
                if m == 1;           
                    nk=size(data1,1);      
                    ij=1;
                    j=2^nk-1;
                    clear c1
                    for k=1:j
                        c= dec2bin(k);
                        a=size(c,2);
                        
                        c3(1:nk)=0;
                        ic=0;
                        for m1=a:-1:1
                            
                            c2=str2num(c(m1));
                            ic=ic+1;
                            c3(nk-ic+1)=c2;
                            if c2 > 0;
                                
                                c1(ij,:)=c3;
                                
                            end
                            
                            
                        end 
                        ij=ij+1;
                        
                    end
                    
                    l=find(c1 == 1); c1(l)=-1; l=find(c1 == 0);c1(l)=1;
                end
                
                
                
                data1_perm=data1-data2;      
                
                
                
                
                for kj=1:nk
                    data1_perm(kj,:)=c1(m,kj).*data1_perm(kj,:);
                end
                PA_perm(1:K) = mean(data1_perm);
                std_PA_perm(1:K)= std(data1_perm);            
                
            end 
            
            
            
        else
            % signe de data1_perm  aléatoire
            clear g
            j=randperm(size(data1,1));
            lm=find (j(1:JA) > JA/2); 
            g(1:JA)=1;
            g(lm)=-1;
            
            g=g';
            
            
            for k=1:K
                data1_perm(:,k)=data1(:,k)'.*g;   
            end
            PA_perm(1:K) = mean(data1_perm);
            std_PA_perm(1:K)= std(data1_perm);
            
            
        end
        
        
        
        STD_PA_PB_perm=sqrt(std_PA_perm.*std_PA_perm./JA);
        % T{m} =abs( PA_perm ) ./ (STD_PA_PB_perm);
        % if abs( PA_perm) > 1.e-06;     %  petite valeurs on néglige
       % T{m} =abs( PA_perm ) ./ (STD_PA_PB_perm);
         T{1} =abs( PA_perm ) ./ (STD_PA_PB_perm);
        % else
        % T{m} =abs( PA_perm) ;
        %end
        % This is the statistics you want to test whether it departs from 0
        
    end  
   % Smax(m)=max(T{m}); 
     Smax(m)=max(T{1}); 
    
    
 %   ll=find(T{m} >= Ps); 
     ll=find(T{1} >= Ps);     Nb_H0(ll)=Nb_H0(ll)+1;
 
end
delete(hwait) 
%Test without p-values-----------------------------------------------------------
S_nop_sort = sort(Smax);
thersh_FWER(1:4)=0;
thersh_FWER(1) = S_nop_sort(floor((length(S_nop_sort)*(1-pthres(1)))));
thersh_FWER(2) = S_nop_sort(floor((length(S_nop_sort)*(1-pthres(2)))));
thersh_FWER(3) = S_nop_sort(floor((length(S_nop_sort)*(1-pthres(3)))));
thersh_FWER(4) = S_nop_sort(floor((length(S_nop_sort)*(1-pthres(4)))));


 for k=1:K; ll=find( Smax >=T_orig (k) ); Pb_max(k)=(size(ll,2)+1)/(M+1); end
Pb_max=sort(Pb_max);

for k=1:K; ll=find( Smax >=T_orig (k) ); Pb_max(k)=(size(ll,2)+1)/(M+1); end
Pb_max=sort(Pb_max);

for k=1:K; Pb_bonf(k)=1-(1-0.05)^k; end
%-Calculate FDR p values
%

%-----------------------------------------------------------------------
Ps=(Nb_H0+1)/(M+1);
p_FDR=zeros(size(Ps));
nd=ndims(Ps);
N(1:nd)=size(Ps);

%psize=prod(N);
Ps=Ps(:);
S=size(Ps,1);
%S=size(Ps,2);
l=find(Ps > 1); Ps(l)=1;
Ps_old=Ps;
[Ps,indice_Ps]=sort(Ps_old);
[PP,ii_ind]=sort(indice_Ps);
cV=1;
% figure
% plot(Ps,'g')
% hold on
% plot(Pb_max,'c')
% plot(Pb_bonf,'b')

%FDR_option== 'BH' ;
Qs  =Ps*S./(1:S)'*cV;   alphat=(1:S)'/(S) ;%-"Corrected" p-values
 thresh_FDR_BH=threshold_FDR(Qs,T_orig,pthres,alphat,ii_ind,nd,N) ;
 %FDR_option== 'BY' ;
 sum_k=0;
 for k=1:S; sum_k=sum_k+1/k ; Qs(k)=Ps(k)*S*sum_k/k;alphat(k)=k/(sum_k*S) ;end
 thresh_FDR_BY=threshold_FDR(Qs,T_orig,pthres,alphat,ii_ind,nd,N) ;
%Sk_FDR=reshape(Ps_cor(ii_ind),N(1:nd));
%Sk_FDR=Qs(ii_ind);
%  reformat the output array
%Sk_FDR=reshape(Sk_FDR,N_data1(2:Siz));

if Siz > 2; T_orig_ori=reshape(T_orig_ori,N_data1(2:Siz));  end
 
%  IN pp
%figure;pcolor(-abs(T_orig_ori));
%colormap(bone); colorbar; shading flat



k1=1;
for k=2:4
    if (thersh_FWER(k) <= thersh_FWER(k-1)); break; end
    k1=k1+1;
end

thersh_FWER=thersh_FWER(1:k1);
pthres=pthres(1:k1);

%apply this threshold to the T_orig statistic   in progress

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   END of the procedure before

return

function thresh_FDR = threshold_FDR(Qs,T_orig,pthres,alphat,ii_ind,nd,N)
    
S=size(Qs,1);
l=find(Qs > 1); Qs(l)=1;
Ps_cor   = zeros(size(Qs));
for i=1:length(Qs)
    % Find "adjusted" p-values
    Ps_cor(i) = min(Qs(i:S));  %-"Adjusted" p-values
end



%plot(alphat)
% plot(Ps_cor,'r')

TT=sort(T_orig);
TT=TT(end:-1:1);
for k=1:size(pthres,2)
   % [a]=find(Qs < pthres(k));
    [a]=find(Ps_cor <= pthres(k));
    if isempty(a) == 0 ;thresh_FDR(k)=TT(max(a)); else; thresh_FDR(k)=0; end
end  

return