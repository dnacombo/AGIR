function varargout = mylimo_random_robust(varargin)

% This function makes the result files for the random effects of various tests
% as well as organizes and makes files for boostrap. Tt is interfaced with
% limo_random_effect which itself interfaces with the user to select and pass
% the data in the appropriate format. Limo_random_robust calls low level
% functions like limo_ttest to perform the actual computation.
%
% FORMAT
%
% varargout{1} outputs the stats
%
% limo_random_robust(1,y,parameter number,nboot)
%                    1 = a one-sample t-test
%                    y = data (dim electrodes, frames, subjects)
%                    parameter number = describe which parameters is analysed (e.g. 1)
%                    nboot = nb of resamples
%
% limo_random_robust(2,y1,y2,parameter number,nboot);
%                    2 = two samples t-test
%                    y1 = data (dim electrodes, frames, subjects)
%                    y2 = data (dim electrodes, frames, subjects)
%                    parameter number = describe which parameters is analysed (e.g. 1)
%                    nboot = nb of resamples
%
% limo_random_robust(3,y1,y2,nboot);
%                    3 = paired t-test
%                    y1 = data (dim electrodes, frames, subjects)
%                    y2 = data (dim electrodes, frames, subjects)
%                    nboot = nb of resamples
%
% limo_random_robust(4,y,X,parameter number,nboot);
%                    4 = regression analysis
%                    y = data (dim electrodes, frames, subjects)
%                    X = design matrix (zscored regressors + cst term)
%                    parameter number = describe which parameters is analysed (e.g. 1)
%                    nboot = nb of resamples
%
% limo_random_robust(5,y,X,gp_nb,cov_nb,nb_items,nboot)
%                    5 = N-way ANOVA/ANCOVA
%                    y = data (dim electrodes, frames, subjects)
%                    nboot = nb of resamples
%
% limo_random_robust(6,y,gp,factor_levels,nboot)
%                    6 = Repeated measures ANOVA/ANCOVA using multivariate approach
%                    y = data (dim electrodes, frames, subjects, measures)
%                    gp = a vector defining gps
%                    factor_levels = a vector specifying the levels of each repeated measure factor
%                    nboot = nb of resamples
%
% limo_random_robust(7,y,gp,factor_levels,nboot)
%                    7 = Repeated measures ANOVA/ANCOVA using univariate approach
%                    y = data (dim electrodes, frames, subjects, measures)
%                    gp = a vector defining gps (only same sample size)
%                    factor_levels = a vector specifying the levels of each repeated measure factor
%                    nboot = nb of resamples
%
% See also LIMO_TTEST LIMO_GLM LIMO_REP_ANOVA

% -----------------------------
%  Copyright (C) LIMO Team 2010

% v1: Cyril Pernet and Guillaume Rousselet 24-08-2009
% v2: Cyril Pernet 12-07-2010
% 25-08-2010: GAR fixed bug in one-sample bootstrap + added new NaN check of boot indices
% 29-08-2010:  Cyril Pernet and Guillaume Rousselet ANOVAs sorted out (boot_index per cell)
% 10-09-2010: Cyril Pernet - made case 6 / 7 so that we either use
%             multivariate stats or univaraite ones for rep \ANOVA

type  = varargin{1};
rand('state',sum(100*clock));
alpha = .05; % used for a basic computation like in t-test but data aren't thresholded here

switch type
    %--------------------------------------------------------------------------
    %            One Sample t-test     //     bootstrap-t method
    %--------------------------------------------------------------------------
    case {1}
        
        data      = varargin{2};
        parameter = varargin{3};
        nboot     = varargin{4};
        clear varargin
        
        % ------------------------------------------------
        % make a one_sample file per parameter (electrodes, frames, [mean value, std, nb_subjects, t, p])
        one_sample = single(NaN(size(data,1), size(data,2), 5));
        %         name = sprintf('one_sample_parameter_%g',parameter);
        s = '';
        for electrode = 1:size(data,1) % run per electrode because we have to remove NaNs
            for i = 1:numel(s)
                fprintf('\b');
            end
            s = sprintf('analyse parameter %g electrode %g',parameter,electrode);
            fprintf(s);
            tmp = squeeze(data(electrode,:,:)); Y = tmp(:,find(~isnan(tmp(1,:))));
            [one_sample(electrode,:,1),dfe,ci,one_sample(electrode,:,2),one_sample(electrode,:,3),one_sample(electrode,:,4),one_sample(electrode,:,5)] = limo_ttest(1,Y,0,alpha);
            clear tmp Y
            if electrode == size(data,1)
                fprintf('\n');
            end
        end
        stats.m = one_sample(:,:,1);
        stats.sd = one_sample(:,:,2);
        stats.n = one_sample(:,:,3);
        stats.t = one_sample(:,:,4);
        stats.p = one_sample(:,:,5);
        varargout{1} = stats;
        varargout{2} = one_sample;
        
        %         save ([name],'one_sample','stats')
        
        % ------------------------------------------------
        % Bootstrap
        if nboot > 0
            
            % create a boot one_sample file to store data under H0 and H1
            boot_one_sample = single(NaN(size(data,1), size(data,2),3,nboot)); % stores T and p values for each boot under H1 then T H0 (last dim)
            boot_name = sprintf('boot_one_sample_parameter_%g',parameter);
            
            % create centered data to estimate H0
            centered_data = data - repmat(nanmean(data,3),[1 1 size(data,3)]);
            
            % create an index to use across all electrodes and frames
            fprintf('\nmaking random table ...\n')
            if size(data,1)==1
                chdata(1,:)=squeeze(data(:,1,:));
            else
                chdata=squeeze(data(:,1,:)); % used to check data for NaNs
            end
            
            B=1;boot_index=zeros(size(data,3),nboot);
            while B~=nboot+1
                tmp = ceil(rand(1,size(data,3)).*size(data,3));
                if length(unique(tmp)) ~= 1 && min(sum(~isnan(chdata(:,tmp)),2)) > 1;
                    boot_index(:,B) = tmp;
                    B=B+1;
                end
            end
            clear chdata
            tmp_boot1 = single(NaN(size(data,2),nboot));
            tmp_boot2 = single(NaN(size(data,2),nboot));
            tmp_boot3 = single(NaN(size(data,2),nboot));
            s = '';
            for electrode = 1:size(data,1)
                
                for B=1:nboot
                    % get results under H1
                    if ~rem(B,100)
                        for i = 1:numel(s)
                            fprintf('\b');
                        end
                        s = sprintf('analyse parameter %g electrode %g boot sample %g',parameter,electrode,B);
                        fprintf(s);
                    end
                    tmp = squeeze(data(electrode,:,boot_index(:,B))); Y = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    [m,dfe,ci,sd,n,tmp_boot1(:,B),p] = limo_ttest(1,Y,0,alpha);
                    % get the result under H0
                    tmp = squeeze(centered_data(electrode,:,boot_index(:,B))); Y2 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    [m,dfe,ci,sd,n,tmp_boot2(:,B),tmp_boot3(:,B)] = limo_ttest(1,Y2,0,alpha);
                    clear Y Y2
                end
                boot_one_sample(electrode,:,1,:) = tmp_boot1; % t values H1
                boot_one_sample(electrode,:,2,:) = tmp_boot2; % t values H0
                boot_one_sample(electrode,:,3,:) = tmp_boot3; % p values H0
                if electrode == size(data,1)
                    fprintf('\n');
                end
                
            end % closes for electrode
            stats.boottH1 = squeeze(boot_one_sample(:,:,1,:));
            stats.boottH0 = squeeze(boot_one_sample(:,:,2,:));
            stats.bootpH0 = squeeze(boot_one_sample(:,:,3,:));
            varargout{1} = stats;
            varargout{3} = boot_one_sample;
            %             save ([boot_name],'boot_one_sample','stats');
        else
            varargout{3} = [];
        end % closes if nboot > 0
        
        
        %--------------------------------------------------------------------------
        %            Two Samples t-test     //     percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {2}
        
        data1      = varargin{2};
        data2      = varargin{3};
        parameter  = varargin{4};
        nboot      = varargin{5};
        clear varargin
        error('don''t forget to deal with varargouts.')
        % ------------------------------------------------
        % make a two_samples file per parameter (electrodes, frames, [mean diff, dfe, std gp1, std gp2, t, p])
        two_samples = single(NaN(size(data1,1), size(data1,2),6));
        name = sprintf('two_samples_ttest_parameter_%g',parameter);
        
        for electrode = 1:size(data1,1) % run per electrode because we have to remove NaNs
            fprintf('analyse parameter %g electrode %g',parameter, electrode); disp(' ');
            tmp = squeeze(data1(electrode,:,:)); Y1 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
            tmp = squeeze(data2(electrode,:,:)); Y2 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
            [two_samples_ttest(electrode,:,1),two_samples_ttest(electrode,:,2), ci, std, n, two_samples_ttest(electrode,:,5), two_samples_ttest(electrode,:,6)] = limo_ttest(2,Y1,Y2,alpha); clear Y1 Y2
            two_samples_ttest(electrode,:,3) = std(:,1)';
            two_samples_ttest(electrode,:,4) = std(:,2)';
        end
        save ([name],'two_samples_ttest')
        
        % ------------------------------------------------
        if nboot > 0
            
            % create a ttest file under H1/H0
            boot_two_samples_ttest = single(NaN(size(data1,1), size(data1,2), 3,nboot)); % stores differences, T and p values for each boot
            boot_name = sprintf('boot_two_samples_ttest_parameter_%g',parameter);
            
            % create an index to use across all electrodes and frames
            disp('making random table ...')
            if size(data1,1)==1
                chdata1(1,:)=squeeze(data1(:,1,:));
                chdata2(1,:)=squeeze(data2(:,1,:));
            else
                chdata1=squeeze(data1(:,1,:)); % used to check data for NaNs
                chdata2=squeeze(data2(:,1,:));
            end
            
            
            B=1;boot_index1=zeros(size(data1,3),nboot);
            while B~=nboot+1
                tmp = ceil(rand(1,size(data1,3)).*size(data1,3));
                if length(unique(tmp)) ~= 1 && min(sum(~isnan(chdata1(:,tmp)),2)) > 1;
                    boot_index1(:,B) = tmp;
                    B=B+1;
                end
            end
            
            B=1;boot_index2=zeros(size(data2,3),nboot);
            while B~=nboot+1
                tmp = ceil(rand(1,size(data2,3)).*size(data2,3));
                if length(unique(tmp)) ~= 1 && min(sum(~isnan(chdata2(:,tmp)),2)) > 1;
                    boot_index2(:,B) = tmp;
                    B=B+1;
                end
            end
            
            tmp_boot1 = single(NaN(size(data1,2),nboot)); tmp_boot2 = single(NaN(size(data2,2),nboot));
            data1_centered = data1 - repmat(nanmean(data1,3),[1 1 size(data1,3)]);
            data2_centered = data2 - repmat(nanmean(data2,3),[1 1 size(data2,3)]);
            
            for electrode = 1:size(data1,1)
                % get results under both H1 and H0
                for B=1:nboot
                    fprintf('analyse parameter %g electrode %g boot sample %g',parameter, electrode,B); disp(' ');
                    % compute the difference under H1 for empirical p values and percentile bootstrap CI
                    tmp = squeeze(data1(electrode,:,boot_index1(:,B))); Y1 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    tmp = squeeze(data2(electrode,:,boot_index2(:,B))); Y2 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    boot_two_samples_ttest(electrode,:,1,B) = nanmean(Y1,2)-nanmean(Y2,2);
                    clear Y1 Y2
                    
                    % to test H0 (cluster) do the two samples t-test on centered data and estimate
                    % when by chance (bootstrap) we have something different from 0
                    tmp = squeeze(data1_centered(electrode,:,boot_index1(:,B))); Y1 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    tmp = squeeze(data2_centered(electrode,:,boot_index2(:,B))); Y2 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    [m,dfe, ci,s,n,boot_two_samples_ttest(electrode,:,2,B),boot_two_samples_ttest(electrode,:,3,B)] = limo_ttest(2,Y1,Y2,alpha);
                    clear Y1 Y2
                end % bootstrap
            end % electrode
            save ([boot_name],'boot_two_samples_ttest');
        end % closes if nboot > 0
        disp('two samples t-test done')
        
        
        %--------------------------------------------------------------------------
        %            Paired t-test     //     percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {3}
        
        data1      = varargin{2};
        data2      = varargin{3};
        nboot      = varargin{4};
        clear varargin
        error('don''t forget to deal with varaergouts.')
        
        % ------------------------------------------------
        % make a paired_ttest file per parameter (electrodes, frames, [mean difference, std, nb_subjects, t, p])
        paired_ttest = single(NaN(size(data1,1), size(data1,2),5));
        
        for electrode = 1:size(data1,1) % run per electrode because we have to remove NaNs
            fprintf('analyse electrode %g',electrode); disp(' ');
            tmp = squeeze(data1(electrode,:,:)); Y1 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
            tmp = squeeze(data2(electrode,:,:)); Y2 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
            [paired_ttest(electrode,:,1),dfe, ci,paired_ttest(electrode,:,2), paired_ttest(electrode,:,3), paired_ttest(electrode,:,4),paired_ttest(electrode,:,5)] = limo_ttest(1,Y1,Y2,alpha); clear Y1 Y2
        end
        save paired_ttest paired_ttest
        
        % ------------------------------------------------
        % get results under H1 and H0
        if nboot > 0
            
            % create a paired ttest under H0
            boot_paired_ttest = single(NaN(size(data1,1), size(data1,2), 3,nboot)); % stores differences, T and p values for each boot
            
            % create an index to use across all electrodes and frames
            disp('making random table ...')
            if size(data1,1)==1
                chdata(1,:)=squeeze(data1(:,1,:));
            else
                chdata=squeeze(data1(:,1,:)); % used to check data for NaNs
            end
            
            B=1;boot_index=zeros(size(data1,3),nboot);
            while B~=nboot+1
                tmp = ceil(rand(1,size(data1,3)).*size(data1,3)); % resample subjects
                if length(unique(tmp)) ~= 1 && min(sum(~isnan(chdata(:,tmp)),2)) > 1;
                    boot_index(:,B) = tmp;
                    B=B+1;
                end
            end
            
            data1_centered = data1 - repmat(nanmean(data1,3),[1 1 size(data1,3)]);
            data2_centered = data2 - repmat(nanmean(data2,3),[1 1 size(data2,3)]);
            
            for electrode = 1:size(data1,1)
                for B=1:nboot
                    fprintf('analyse electrode %g boot sample %g',electrode,B); disp(' ');
                    % compute the difference under H1 for corrected p values and CI
                    tmp = squeeze(data1(electrode,:,boot_index(:,B))); Y1 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    tmp = squeeze(data2(electrode,:,boot_index(:,B))); Y2 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    boot_paired_ttest(electrode,:,1,B) = nanmean((Y1 - Y2),2);
                    % compute the paired t-test under H0 for the cluster analysis
                    tmp = squeeze(data1_centered(electrode,:,boot_index(:,B))); Y1 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    tmp = squeeze(data2_centered(electrode,:,boot_index(:,B))); Y2 = tmp(:,find(~isnan(tmp(1,:)))); clear tmp
                    [m,dfe,ci,s,n,boot_paired_ttest(electrode,:,2,B),boot_paired_ttest(electrode,:,3,B)] = limo_ttest(1,Y1,Y2,alpha);
                    clear Y1 Y2
                end
            end
            save boot_paired_ttest boot_paired_ttest
        end
        disp('paired t-test done')
        
        
        
        %--------------------------------------------------------------------------
        %                    Regression  // percentile bootstrap under H0 and H1
        %--------------------------------------------------------------------------
    case {4}
        error('don''t forget to deal with varaergouts.')
        
        data       = varargin{2};
        regressors = varargin{3}; % the predictors across subjects like e.g. age
        parameter  = varargin{4}; % the parameters from 1st level matrices the regression is computed on, coded as e.g. 1
        nboot      = varargin{5};
        clear varargin
        
        % ------------------------------------------------
        % make the Betas and regression (results) files per parameter (electrodes, frames, nb_regressors, [F / p])
        Betas = single(NaN(size(data,1), size(data,2), size(regressors,2)));
        regression = single(NaN(size(data,1), size(data,2), size(regressors,2)-1, 2));
        name = sprintf('regression_parameter_%g',parameter);
        Res = single(NaN(size(data)));
        
        % compute the regression
        for electrode = 1:size(data,1)
            fprintf('analyse parameter %g electrode %g',parameter, electrode); disp(' ');
            
            tmp = squeeze(data(electrode,:,:));
            Y = tmp(:,find(~isnan(tmp(1,:))))';
            X = regressors(find(~isnan(tmp(1,:)))',:); % Y and X with NaN from Y data removed
            if sum(mean(X)) - 1 ~= 0  % zscore the covariates if needed - note it might be diff. for different electrodes because of the NaNs in Y removed above in Y and X
                X(:,1:end-1) = zscore(X(:,1:end-1));
            end
            
            model = limo_glm(Y,X,0,size(X,2)-1,0,1);
            Betas(electrode,:,:) = model.betas';
            fitted_data = X*model.betas;
            regression(electrode,:,:,1) = model.univariate.continuous.F'; % note the intercept has no F and p values (not coded)
            regression(electrode,:,:,2) = model.univariate.continuous.p';
            good_subjects = find(~isnan(tmp(1,:)));
            Res(electrode,:,good_subjects)  = (Y - fitted_data)';
            clear tmp Y X model fitted_data
        end
        save Betas Betas; save Res Res; save ([name],'regression')
        
        % ----------------------------------------------------------------
        if nboot > 0
            
            % create a regression file under H1/H0
            boot_H1_Betas = single(NaN(size(data,1), size(data,2), size(regressors,2), nboot)); % only stores betas for H1
            boot_H0_Betas = single(NaN(size(data,1), size(data,2), size(regressors,2), nboot));
            boot_H0_regression = single(NaN(size(data,1), size(data,2), size(regressors,2)-1, 2, nboot)); % stores F and p values for each boot
            boot_name_H0 = sprintf('boot_H0_regression_parameter_%g',parameter);
            
            % create an index to use across all electrodes and frames
            disp('making random table ...')
            if size(data,1)==1
                chdata(1,:)=squeeze(data(:,1,:));
            else
                chdata=squeeze(data(:,1,:)); % used to check data for NaNs
            end
            
            B=1;boot_index=zeros(size(data,3),nboot);
            while B~=nboot+1
                tmp = ceil(rand(1,size(data,3)).*size(data,3));
                if length(unique(tmp)) ~= 1 && min(sum(~isnan(chdata(:,tmp)),2)) > 1;
                    boot_index(:,B) = tmp;
                    B=B+1;
                end
            end
            
            for B=1:nboot
                for electrode = 1:size(data,1)
                    fprintf('analyse parameter %g boot sample %g electrode %g',parameter,B,electrode); disp(' ');
                    % compute the effects under H0 for MCC
                    tmp = squeeze(data(electrode,:,boot_index(:,B)));
                    Y = tmp(:,find(~isnan(tmp(1,:))))';
                    X = regressors(find(~isnan(tmp(1,:)))',:); % removes the same rows as Y for size issues but Y has been resampled not X, i.e. Y and X do not match = H0
                    if sum(mean(X)) - 1 ~= 0 % zscoring
                        X(:,1:end-1) = zscore(X(:,1:end-1));
                    end
                    
                    model = limo_glm(Y,X,0,size(X,2)-1,0,1);
                    boot_H0_Betas(electrode,:,:,B) = model.betas';
                    boot_H0_regression(electrode,:,:,1,B) = model.univariate.continuous.F';
                    boot_H0_regression(electrode,:,:,2,B) = model.univariate.continuous.p';
                    
                    % compute the effects under H1 for empirical p values and CI (Wilcox p418)
                    tmp2 = regressors(boot_index(:,B),:); % sample X as Y
                    X2 = tmp2(find(~isnan(tmp(1,:)))',:); % keeps the regressors matched to Y
                    if sum(mean(X2)) - 1 ~= 0 % zscoring
                        X2(:,1:end-1) = zscore(X2(:,1:end-1));
                    end
                    
                    boot_H1_Betas(electrode,:,:,B) = (pinv(X2)*Y)';   % simply use the beta distrib. to test if significant under H1
                    clear tmp Y X model tmp2 X2
                end
            end
            save boot_index boot_index
            save boot_H0_Betas boot_H0_Betas
            save boot_H1_Betas boot_H1_Betas
            save ([boot_name_H0],'boot_H0_regression');
        end
        
        %--------------------------------------------------------------------------
        %                    One-way ANOVA / ANCOVA
        %--------------------------------------------------------------------------
    case {5}
        error('don''t forget to deal with varaergouts.')
        
        data          = varargin{2};
        design_matrix = varargin{3};
        nb_gp         = varargin{4};
        nb_cont       = varargin{5};
        nb_item       = varargin{6};
        nboot         = varargin{7};
        clear varargin
        
        
        %--------------------------------------------------------------------------
        %                    One-way ANOVA - bootstrap by centering data
        %--------------------------------------------------------------------------
        if nb_cont == 0
            
            % make a file per effect (electrodes, frames, [F / p])
            Res = single(NaN(size(data)));
            Betas = single(NaN(size(data,1), size(data,2), nb_gp+1));
            N_way_variance_analysis = single(NaN(size(data,1), size(data,2), 2)); % stores F & p values
            name = sprintf('One_way_variance_analysis');
            
            % compute the ANOVA
            for electrode = 1:size(data,1)
                fprintf('analyse electrode %g', electrode); disp(' ');
                tmp = squeeze(data(electrode,:,:));
                Y = tmp(:,find(~isnan(tmp(1,:))))';
                X = design_matrix(find(~isnan(tmp(1,:)))',:); % removes the same rows as Y
                model = limo_glm(Y,X,nb_gp,nb_cont,sum(X(:,1:nb_gp),1),1);
                Betas(electrode,:,:) = model.betas';
                N_way_variance_analysis(electrode,:,1) = model.univariate.conditions.F';
                N_way_variance_analysis(electrode,:,2) = model.univariate.conditions.p';
                fitted_data = X*model.betas;
                good_subjects = find(~isnan(tmp(1,:)));
                Res(electrode,:,good_subjects)  = (Y - fitted_data)';
                clear tmp Y X model
            end
            save Res Res
            save Betas Betas
            save ([name],'N_way_variance_analysis')
            clear Betas N_way_variance_analysis
            
            
            % ----------------------------------------------------------------
            if nboot > 0
                disp('starting bootstrap ...')
                
                % create an index to use across all electrodes and frames
                disp('making random table...')
                index1=1;index2=nb_item(1);
                boot_index=zeros(size(data,3),nboot);
                for cel=1:length(nb_item)
                    
                    if size(data,1)==1
                        chdata(1,:)=squeeze(data(1,1,index1:index2));
                    else
                        chdata=squeeze(data(:,1,index1:index2)); % used to check data for NaNs
                    end
                    
                    B=1;
                    while B~=nboot+1
                        tmp = ceil(rand(1,nb_item(cel)).*nb_item(cel)+index1-1);
                        if length(unique(tmp)) ~= 1 && min(sum(~isnan(chdata(:,tmp-index1+1)),2)) > 1;
                            boot_index(index1:index2,B) = tmp'; % works because subjects are stacked by group on dim(3) of data
                            B=B+1;
                        end
                    end
                    
                    if cel<length(nb_item)
                        index1 = index1+nb_item(cel);
                        index2 = index2+nb_item(cel+1);
                    end
                    clear chdata
                end
                
                
                % create a file to store bootstrap Beta values under H1
                boot_H1_Betas = single(NaN(size(data,1), size(data,2), nb_gp+1, nboot));
                for B=1:nboot
                    fprintf('create robust paramters ... boot %g',B); disp(' ');
                    for electrode = 1:size(data,1)
                        tmp = squeeze(data(electrode,:,boot_index(:,B))); % sample from all subjects in each group
                        Y = tmp(:,find(~isnan(tmp(1,:))))'; % remove NaNs
                        X = design_matrix(boot_index(:,B),:); % resample X as Y
                        X = X(find(~isnan(tmp(1,:)))',:); % now removes the same rows as Y
                        boot_H1_Betas(electrode,:,:,B) = (pinv(X)*Y)';
                        clear tmp Y X
                    end
                end
                save boot_H1_Betas boot_H1_Betas ; clear boot_H1_Betas
                
                % create files to store F and p bootstrap under H0
                boot_H0_Betas = single(NaN(size(data,1), size(data,2), nb_gp+1, nboot));
                boot_H0_N_way_variance_analysis = single(NaN(size(data,1), size(data,2), 2, nboot)); % stores F and p values for each boot
                boot_name = sprintf('boot_H0_One_way_variance_analysis');
                
                % center each cell to compute under H0
                centered_data = single(NaN(size(data,1),size(data,2),size(data,3)));
                for cel=1:length(nb_item)
                    index = find(design_matrix(:,cel));
                    centered_data(:,:,index) = data(:,:,index) - repmat(nanmean(data(:,:,index),3),[1 1 size(data(:,:,index),3)]);
                end
                
                % clear memory
                clear data index
                
                % do the bootstrap on centered data
                for B=1:nboot
                    for electrode = 1:size(centered_data,1)
                        fprintf('analyse boot sample %g electrode %g',B,electrode); disp(' ');
                        tmp = squeeze(centered_data(electrode,:,boot_index(:,B))); % sample randomly from each group of subjects
                        Y = tmp(:,find(~isnan(tmp(1,:))))'; % remove NaNs
                        X = design_matrix(boot_index(:,B),:); % resample X as Y
                        X = X(find(~isnan(tmp(1,:)))',:); % now removes the same rows as Y
                        model = limo_glm(Y,X,nb_gp,nb_cont,sum(X(:,1:nb_gp),1),1);
                        boot_H0_Betas(electrode,:,:,B) = model.betas';
                        boot_H0_N_way_variance_analysis(electrode,:,1,B) = model.univariate.conditions.F';
                        boot_H0_N_way_variance_analysis(electrode,:,2,B) = model.univariate.conditions.p';
                        clear tmp Y X model
                    end
                end
                
                % save results
                disp('saving results .. ')
                save ([boot_name],'boot_H0_N_way_variance_analysis');
                save boot_H0_Betas boot_H0_Betas
                
                % save centered_data and boot_index to be used for contrasts
                save centered_data centered_data
                save boot_index boot_index
            end
            disp('ANOVA done')
            
            
            %--------------------------------------------------------------------------
            %                    One-way ANCOVA - bootstrap by breaking the link Y/X
            %--------------------------------------------------------------------------
        else
            
            % make a file per effect (electrodes, frames, nb_regressors, [F / p])
            Betas = single(NaN(size(data,1), size(data,2), nb_gp+nb_cont+1));
            N_way_covariance_analysis = single(NaN(size(data,1), size(data,2), 1+nb_cont, 2));
            name = sprintf('One_way_covariance_analysis');
            
            % compute the ANCOVA
            for electrode = 1:size(data,1)
                fprintf('analyse electrode %g', electrode); disp(' ');
                tmp = squeeze(data(electrode,:,:));
                Y = tmp(:,find(~isnan(tmp(1,:))))';
                X = design_matrix(find(~isnan(tmp(1,:)))',:); % removes the same rows as Y
                if sum(mean(X(:,nb_gp+1:end-1),1)) ~= 0  % zscore the covariates
                    X(:,nb_gp+1:end-1) = zscore(X(:,nb_gp+1:end-1));
                end
                
                model = limo_glm(Y,X,nb_gp,nb_cont,sum(X(:,1:nb_gp),1),1);
                Betas(electrode,:,:) = model.betas';
                N_way_covariance_analysis(electrode,:,1,1) = model.univariate.conditions.F';
                N_way_covariance_analysis(electrode,:,1,2) = model.univariate.conditions.p';
                N_way_covariance_analysis(electrode,:,2:end,1) = model.univariate.continuous.F';
                N_way_covariance_analysis(electrode,:,2:end,2) = model.univariate.continuous.p';
                fitted_data = X*model.betas;
                good_subjects = find(~isnan(tmp(1,:)));
                Res(electrode,:,good_subjects)  = (Y - fitted_data)';
            end
            save Res Res
            save Betas Betas
            save ([name],'N_way_covariance_analysis')
            clear Betas N_way_covariance_analysis
            
            % ----------------------------------------------------------------
            if nboot > 0
                
                % create an index to use across all electrodes and frames
                disp('making random table...')
                index1=1;index2=nb_item(1);
                boot_index=zeros(size(data,3),nboot);
                for cel=1:length(nb_item)
                    
                    if size(data,1)==1
                        chdata(1,:)=squeeze(data(1,1,index1:index2));
                    else
                        chdata=squeeze(data(:,1,index1:index2)); % used to check data for NaNs
                    end
                    
                    B=1;
                    while B~=nboot+1
                        tmp = ceil(rand(1,nb_item(cel)).*nb_item(cel)+index1-1);
                        if length(unique(tmp)) ~= 1 && min(sum(~isnan(chdata(:,tmp-index1+1)),2)) > 1;
                            boot_index(index1:index2,B) = tmp; % works because subjects are stacked by group on dim(3) of data
                            B=B+1;
                        end
                    end
                    
                    if cel<length(nb_item)
                        index1 = index1+nb_item(cel);
                        index2 = index2+nb_item(cel+1);
                    end
                    clear chdata
                end
                
                % create files to store bootstrap under H1
                % -----------------------------------------
                boot_H1_Betas = single(NaN(size(data,1), size(data,2), nb_gp+nb_cont+1, nboot));
                
                for B=1:nboot
                    fprintf('computing paramtewrs under H1 .. boot sample %g %g',B,electrode); disp(' ');
                    for electrode = 1:size(data,1)
                        tmp = squeeze(data(electrode,:,boot_index(:,B))); % sample randomly from all Y
                        Y = tmp(:,find(~isnan(tmp(1,:))))'; % remove NaNs
                        X = design_matrix(boot_index(:,B),:); % resample X as Y
                        X = X(find(~isnan(tmp(1,:)))',:); % removes the same rows as Y (Y and X match)
                        if sum(mean(X(:,nb_gp+1:end-1),1)) ~= 0  % zscore the covariates
                            X(:,nb_gp+1:end-1) = zscore(X(:,nb_gp+1:end-1));
                        end
                        boot_H1_Betas(electrode,:,:,B) = (pinv(X)*Y)';
                        clear tmp Y X2
                    end
                end
                
                save boot_H1_Betas boot_H1_Betas
                clear boot_H1_Beatas
                
                % create files to store bootstrap under H0
                % -----------------------------------------
                boot_H0_Betas = single(NaN(size(data,1), size(data,2), nb_gp+nb_cont+1, nboot));
                boot_H0_N_way_covariance_analysis = single(NaN(size(data,1), size(data,2), 1+nb_cont, 2, nboot)); % stores F and p values for each boot
                boot_name = sprintf('boot_H0_One_way_covariance_analysis');
                
                for B=1:nboot
                    for electrode = 1:size(data,1)
                        fprintf('analyse boot sample %g electrode %g',B,electrode); disp(' ');
                        
                        % do the analysis under H0 (note here H0 is
                        % estimated by breaking the link bertween Y and X
                        % rather that working on centered data)
                        tmp = squeeze(data(electrode,:,boot_index(:,B))); % sample randomly from all Y
                        Y = tmp(:,find(~isnan(tmp(1,:))))'; % remove NaNs
                        X = design_matrix(find(~isnan(tmp(1,:)))',:); % removes the same rows as Y for size issue but X is not resampled, ie Y and X don't match anymore
                        if sum(mean(X(:,nb_gp+1:end-1),1)) ~= 0  % zscore the covariates
                            X(:,nb_gp+1:end-1) = zscore(X(:,nb_gp+1:end-1));
                        end
                        
                        model = limo_glm(Y,X,nb_gp,nb_cont,sum(X(:,1:nb_gp),1),1);
                        boot_H0_Betas(electrode,:,:,B) = model.betas';
                        boot_H0_N_way_covariance_analysis(electrode,:,1,1,B) = model.univariate.conditions.F';
                        boot_H0_N_way_covariance_analysis(electrode,:,1,2,B) = model.univariate.conditions.p';
                        boot_H0_N_way_covariance_analysis(electrode,:,2:end,1,B) = model.univariate.continuous.F';
                        boot_H0_N_way_covariance_analysis(electrode,:,2:end,2,B) = model.univariate.continuous.p';
                        clear X model
                    end
                end
                
                % spend time to load/save but necessary for memory issue
                save boot_index boot_index
                save boot_H0_Betas boot_H0_Betas
                save ([boot_name],'boot_H0_N_way_covariance_analysis');
                clear boot_H1_Betas boot_H0_Betas boot_H0_N_way_covariance_analysis
                
            end
        end
        disp('ANCOVA done')
        
        
        %----------------------------------------------------------------------------------------------
        %                    Repeated Measure ANOVA (multivariate approach) - bootstrap centering data
        %----------------------------------------------------------------------------------------------
    case {6}
        error('don''t forget to deal with varaergouts.')
        
        data          = varargin{2};
        gp_vector     = varargin{3};
        factor_levels = varargin{4};
        nboot         = varargin{5};
        clear varargin
        
        % from the input we know which case to handle
        if unique(gp_vector) == 1
            % one sample
            if length(factor_levels) ==1
                type = 1;
            elseif length(factor_levels) >1
                type = 2;
            end
        else
            % k samples
            if length(factor_levels) ==1
                type = 3;
            elseif length(factor_levels) >1
                type = 4;
            end
        end
        
        
        % make files to be stored
        if type == 1 % one factor
            C                     = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            Rep_ANOVA             = single(NaN(size(data,1),size(data,2),3)); % store mean difference, F and p
            load LIMO;
            LIMO.design.Name{1} = 'Main effect';
            LIMO.design.C = C;
            LIMO.design.averages = factor_levels;
            x = kron(eye(prod(factor_levels)),ones(size(data,3),1));
            LIMO.design.X = [x ones(size(x,1),1)];
            
        elseif type == 2 % many factor
            C = limo_OrthogContrasts(factor_levels);
            for i=1:length(C); m(i) = size(C{i},1); end
            Rep_ANOVA             = single(NaN(size(data,1),size(data,2),length(C),3)); % store mean difference, F and p for each within factor and interactions
            load LIMO; LIMO.design.C = C; LIMO.design.averages = m; index = length(factor_levels)+1;
            for i= 1:length(factor_levels); LIMO.design.name{i} = ['Main effect ' num2str(i)]; end
            for i= 2:length(factor_levels); n = nchoosek([1:length(factor_levels)],i);
                for j=1:size(n,1); LIMO.design.name{index} = ['Interaction ' num2str(n(j,:))]; index = index+1; end; end
            x = kron(eye(prod(factor_levels)),ones(size(data,3),1));
            LIMO.design.X = [x ones(size(x,1),1)];
            
        elseif type == 3 % one factor within and one factor between
            gp_values = unique(gp_vector); k = length(gp_values); X = single(NaN(size(gp_vector,1),k+1));
            for g =1:k; X(:,g) = gp_vector == gp_values(g); end; X(:,end) = 1; % design matrix for gp effects
            G = zeros(size(gp_vector,1)*prod(factor_levels),k+1); G_vector = repmat(gp_vector,prod(factor_levels),1);
            index1 = 1; index2 = max(find(X(:,1))*prod(factor_levels));
            for g =1:k; G(index1:index2,g) = 1; index1 = index2+1; index2 = max(find(X(:,g+1))*prod(factor_levels)); end; G(:,end) = 1; % G is form the display only
            C                             = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            Rep_ANOVA                     = single(NaN(size(data,1),size(data,2),3));
            Rep_ANOVA_Gp_effect           = single(NaN(size(data,1),size(data,2),3));
            Rep_ANOVA_Interaction_with_gp = single(NaN(size(data,1),size(data,2),3));
            load LIMO; LIMO.design.C = C; LIMO.design.Name = 'Main effect';
            x = []; for g=1:k; x = [x;kron(eye(prod(factor_levels)),ones(sum(X(:,g)),1))]; end
            LIMO.design.X = [G(:,1:end-1) x G(:,end)];
            
        elseif type == 4 % many factors within and one factor between
            gp_values = unique(gp_vector); k = length(gp_values); X = single(NaN(size(gp_vector,1),k+1));
            for g =1:k; X(:,g) = gp_vector == gp_values(g); end; X(:,end) = 1; % design matrix for gp effects
            G = zeros(size(gp_vector,1)*prod(factor_levels),k+1); G_vector = repmat(gp_vector,prod(factor_levels),1);
            index1 = 1; index2 = max(find(X(:,1))*prod(factor_levels));
            for g =1:k; G(index1:index2,g) = 1; index1 = index2+1; index2 = max(find(X(:,g+1))*prod(factor_levels)); end; G(:,end) = 1; % G is form the display only
            C = limo_OrthogContrasts(factor_levels);
            Rep_ANOVA                     = single(NaN(size(data,1),size(data,2),length(C),3));
            Rep_ANOVA_Gp_effect           = single(NaN(size(data,1),size(data,2),3));
            Rep_ANOVA_Interaction_with_gp = single(NaN(size(data,1),size(data,2),length(C),3));
            load LIMO; LIMO.design.C = C; for i= 1:length(factor_levels); LIMO.design.name{i} = ['Main effect ' num2str(i)]; end
            index = length(factor_levels)+1;
            for i= 2:length(factor_levels);
                n = nchoosek([1:length(factor_levels)],i);
                for j=1:size(n,1)
                    LIMO.design.name{index} = ['Interaction ' num2str(n(j,:))];
                    index = index+1;
                end
            end
            x = []; for g=1:k; x = [x;kron(eye(prod(factor_levels)),ones(sum(X(:,g)),1))]; end
            LIMO.design.X = [G(:,end-1) x G(:,end)];
        end
        
        % check the design with user
        figure('Name','Design matrix'); set(gcf,'Color','w'); imagesc(LIMO.design.X);
        colormap('gray');  title('ANOVA model','FontSize',16);xlabel('regressors');
        ylabel('subjects'); drawnow;
        go = questdlg('start the analysis?');
        switch go
            case {'No' 'Cancel'}
                return
        end
        save LIMO LIMO; clear LIMO
        
        % do the analysis
        for electrode = 1:size(data,1)
            fprintf('analyse electrode %g ...', electrode); disp(' ');
            for frame = 1:size(data,2)
                tmp = squeeze(data(electrode,frame,:,:));
                Y = tmp(find(~isnan(tmp(:,1))),:);
                gp = gp_vector(find(~isnan(tmp(:,1))),:);
                if type == 3 || type == 4
                    XB = X(find(~isnan(tmp(:,1))),:);
                end
                
                if type == 1
                    result = limo_rep_anova(Y,gp,factor_levels,C,cov(Y));
                    Rep_ANOVA(electrode,frame,1)                       = nanmean(C*Y');
                    Rep_ANOVA(electrode,frame,2:3)                     = [result.F;result.p];
                    
                elseif type == 2
                    result = limo_rep_anova(Y,gp,factor_levels,C,cov(Y));
                    Rep_ANOVA(electrode,frame,:,2)                     = result.F;
                    Rep_ANOVA(electrode,frame,:,3)                     = result.p;
                    for effect = 1:size(C,2)
                        c           = C{effect};
                        avg(effect) = nanmean(c*Y');
                    end
                    Rep_ANOVA(electrode,frame,:,1) = avg;
                    
                elseif type == 3
                    result = limo_rep_anova(Y,gp,factor_levels,C,XB);
                    Rep_ANOVA(electrode,frame,2:3)                       = [result.repeated_measure.F;result.repeated_measure.p];
                    Rep_ANOVA_Gp_effect(electrode,frame,2:3)             = [result.gp.F;result.gp.p];
                    Rep_ANOVA_Interaction_with_gp(electrode,frame,2:3)   = [result.interaction.F;result.interaction.p];
                    yp  = nanmean(Y,1)'; Rep_ANOVA(electrode,frame,1) = nanmean(C*yp);
                    for g=1:k
                        yg(g) = nanmean(nanmean(Y(XB(:,g)==1,:))); % mean per gp
                    end
                    Rep_ANOVA_Gp_effect(electrode,frame,1) = nanmean([eye(length(yg)-1) ones(length(yg)-1,1).*-1]*yg');
                    I = (C*Y')'; Rep_ANOVA_Interaction_with_gp(electrode,frame,1) = nanmean(nanmean(I));
                    
                elseif type == 4
                    result = limo_rep_anova(Y,gp,factor_levels,C,X);
                    Rep_ANOVA(electrode,frame,:,1)                     = result.repeated_measure.means;
                    Rep_ANOVA(electrode,frame,:,2)                     = result.repeated_measure.F;
                    Rep_ANOVA(electrode,frame,:,3)                     = result.repeated_measure.p;
                    Rep_ANOVA_Gp_effect(electrode,frame,:)             = [result.gp.mean;result.gp.F;result.gp.p];
                    Rep_ANOVA_Interaction_with_gp(electrode,frame,:,1) = rresult.interaction.means;
                    Rep_ANOVA_Interaction_with_gp(electrode,frame,:,2) = result.interaction.F;
                    Rep_ANOVA_Interaction_with_gp(electrode,frame,:,3) = result.interaction.p;
                    for effect = 1:length(C)
                        c   = C{effect};  I = (c*Y')';
                        avg(effect) = nanmean(nanmean(c*I));
                        avg2(effect) = nanmean(nanmean(I));
                    end
                    Rep_ANOVA(electrode,frame,:,1) = avg;
                    Rep_ANOVA_Interaction_with_gp(electrode,frame,:,1) = avg2;
                    for g=1:k
                        yg(g) = nanmean(nanmean(Y(XB(:,g)==1,:))); % mean per gp
                    end
                    Rep_ANOVA_Gp_effect(electrode,frame,1) = nanmean([eye(length(yg)-1) ones(length(yg)-1,1).*-1]*yg');
                    
                end
                clear tmp Y gp result
            end
        end
        
        % save stuff
        if type == 1 || type ==2
            save Rep_ANOVA Rep_ANOVA; clear Rep_ANOVA
        else
            save Rep_ANOVA Rep_ANOVA; clear Rep_ANOVA
            save Rep_ANOVA_Gp_effect Rep_ANOVA_Gp_effect; clear Rep_ANOVA_Gp_effect
            save Rep_ANOVA_Interaction_with_gp Rep_ANOVA_Interaction_with_gp; clear Rep_ANOVA_Interaction_with_gp
        end
        
        
        % ----------------------------------------------------------------
        if nboot > 0
            
            % create files to store bootstrap under H1 and H0
            disp('making bootstrap files ...')
            if type ==1
                boot_H0_Rep_ANOVA                     = single(NaN(size(data,1),size(data,2),3,nboot));
            elseif type == 2
                boot_H0_Rep_ANOVA                     = single(NaN(size(data,1),size(data,2),length(C),3,nboot));
            elseif type == 3
                boot_H0_Rep_ANOVA                     = single(NaN(size(data,1),size(data,2),3));
                boot_H0_Rep_ANOVA_Gp_effect           = single(NaN(size(data,1),size(data,2),3));
                boot_H0_Rep_ANOVA_Interaction_with_gp = single(NaN(size(data,1),size(data,2),3));
            else
                boot_H0_Rep_ANOVA                     = single(NaN(size(data,1),size(data,2),length(C),3));
                boot_H0_Rep_ANOVA_Gp_effect           = single(NaN(size(data,1),size(data,2),3));
                boot_H0_Rep_ANOVA_Interaction_with_gp = single(NaN(size(data,1),size(data,2),length(C),3));
            end
            
            % the data have to be centered (H0) for each cell
            centered_data = single(NaN(size(data,1),size(data,2),size(data,3),size(data,4)));
            
            if type ==3 || type ==4
                ss=sum(X(:,1:end-1));
            else
                ss = size(centered_data,3);
            end
            
            gp_index1 = 1;
            nb_conditions = prod(factor_levels);
            for gp=1:length(ss)
                gp_index2 = gp_index1-1+ss(gp);
                for condition=1:nb_conditions
                    avg = repmat(nanmean(data(:,:,gp_index1:gp_index2,condition),3),[1 1 length(gp_index1:gp_index2)]);
                    centered_data(:,:,gp_index1:gp_index2,condition) = data(:,:,gp_index1:gp_index2,condition) - avg; clear avg
                end
                gp_index1 = gp_index1+ss(gp);
            end
            clear gp
            save centered_data centered_data
            
            
            % create an index to use across all electrodes and frames
            % (different per gp but identical across conditions)
            disp('making random table...')
            index1=1;index2=ss(1);
            boot_index=single(NaN(size(centered_data,3),nboot));
            for cel=1:length(ss)
                if size(centered_data,1)==1
                    chdata(1,:)=squeeze(centered_data(1,1,index1:index2,1)); % used to check data for NaNs - only take 1 condition since NaNs are on electrodes
                else
                    chdata=squeeze(centered_data(:,1,index1:index2,1));
                end
                
                B=1;
                while B~=nboot+1
                    tmp = (ceil(rand(1,ss(cel)).*ss(cel)+index1-1));
                    if length(unique(tmp)) ~= 1 && min(sum(~isnan(chdata(:,tmp-index1+1)),2)) > 1;
                        boot_index(index1:index2,B) = tmp; % works because subjects are stacked by group on dim(3) of data
                        B=B+1;
                    end
                end
                
                if cel<length(ss)
                    index1 = index1+ss(cel);
                    index2 = index2+ss(cel+1);
                end
                clear chdata
            end
            
            
            % compute bootstrap under H1 for means and under H0 for F and p
            for B=1:nboot
                for electrode = 1:size(centered_data,1)
                    fprintf('bootstrap %g electrode %g ...', B, electrode); disp(' ');
                    % H1
                    tmp = squeeze(data(electrode,:,boot_index(:,B),:));
                    Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
                    gp = gp_vector(find(~isnan(tmp(1,:,1))));
                    if type == 3 || type == 4
                        XB = X(find(~isnan(tmp(1,:,1))));
                    end
                    
                    for frame = 1:size(centered_data,2)
                        y = squeeze(Y(frame,:,:));
                        if type == 1
                            boot_H0_Rep_ANOVA(electrode,frame,1,B) = nanmean(C*y');
                        elseif type == 2
                            for effect = 1:size(C,2)
                                c           = C{effect};
                                avg(effect) = nanmean(c*y');
                            end
                            boot_H0_Rep_ANOVA(electrode,frame,:,1,B) = avg;
                        elseif type == 3
                            yp  = nanmean(y,1)';
                            boot_H0_Rep_ANOVA(electrode,frame,1,B) = nanmean(C*yp);
                            for g=1:k
                                yg = nanmean(nanmean(y(XB(:,g)==1,:))); % mean per gp
                            end
                            boot_H0_Rep_ANOVA_Gp_effect(electrode,frame,1,B) = yg;
                            I = (C*y')';
                            boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,frame,1,B) = nanmean(nanmean(I));
                        elseif type == 4
                            for effect = 1:length(C)
                                c   = C{effect};  I = (c*y')';
                                avg(effect) = nanmean(nanmean(c*y));
                                avg2(effect) = nanmean(nanmean(I));
                            end
                            boot_H0_Rep_ANOVA(electrode,frame,:,1,B) = avg;
                            boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,frame,:,1,B) = avg2;
                            for g=1:k
                                yg = nanmean(nanmean(y(XB(:,g)==1,:))); % mean per gp
                            end
                            boot_H0_Rep_ANOVA_Gp_effect(electrode,frame,1,B) = yg;
                        end
                        clear y
                    end
                    clear XB Y gp tmp
                    
                    % H0
                    tmp = squeeze(centered_data(electrode,:,boot_index(:,B),:));
                    Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
                    gp = gp_vector(find(~isnan(tmp(1,:,1))));
                    if type == 3 || type == 4
                        XB = X(find(~isnan(tmp(1,:,1))));
                    end
                    
                    for frame = 1:size(centered_data,2)
                        y = squeeze(Y(frame,:,:));
                        if type == 1
                            result = limo_rep_anova(y,gp,factor_levels,C,cov(y));
                            boot_H0_Rep_ANOVA(electrode,frame,2:3,B)                     = [result.F;result.p];
                        elseif type == 2
                            result = limo_rep_anova(y,gp,factor_levels,C,cov(y));
                            boot_H0_Rep_ANOVA(electrode,frame,:,2,B)                     = result.F;
                            boot_H0_Rep_ANOVA(electrode,frame,:,3,B)                     = result.p;
                        elseif type == 3
                            result = limo_rep_anova(y,gp,factor_levels,C,XB);
                            boot_H0_Rep_ANOVA(electrode,frame,2:3,B)                     = [result.repeated_measure.F;result.repeated_measure.p];
                            boot_H0_Rep_ANOVA_Gp_effect(electrode,frame,2:3,B)           = [result.gp.F;result.gp.p];
                            boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,frame,2:3,B) = [result.interaction.F;result.interaction.p];
                        elseif type == 4
                            result = limo_rep_anova(y,gp,factor_levels,C,XB);
                            boot_H0_Rep_ANOVA(electrode,frame,:,2,B)                     = result.repeated_measure.F;
                            boot_H0_Rep_ANOVA(electrode,frame,:,3,B)                     = result.repeated_measure.p;
                            boot_H0_Rep_ANOVA_Gp_effect(electrode,frame,:,B)             = [result.gp.F;result.gp.p];
                            boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,frame,:,2,B) = result.interaction.F;
                            boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,frame,:,3,B) = result.interaction.p;
                        end
                        clear y result
                    end
                    clear XB Y gp tmp
                end
            end
            
            % save boot_H1_Betas boot_H1_Betas; clear boot_H1_Betas
            if type == 1 || type ==2
                save boot_H0_Rep_ANOVA boot_H0_Rep_ANOVA; clear boot_H0_Rep_ANOVA
            else
                save boot_H0_Rep_ANOVA boot_H0_Rep_ANOVA; clear boot_H0_Rep_ANOVA
                save boot_H0_Rep_ANOVA_Gp_effect boot_H0_Rep_ANOVA_Gp_effect; clear boot_H0_Rep_ANOVA_Gp_effect
                save boot_H0_Rep_ANOVA_Interaction_with_gp boot_H0_Rep_ANOVA_Interaction_with_gp; clear boot_H0_Rep_ANOVA_Interaction_with_gp
            end
        end
        
        
        %---------------------------------------------------------------------------------------------
        %                    Repeated Measure ANOVA (univariate apoproach) - bootstrap centering data
        %----------------------------------------------------------------------------------------------
        
    case {7}
        error('don''t forget to deal with varaergouts.')
        
        data          = varargin{2};
        gp_vector     = varargin{3};
        factor_levels = varargin{4};
        nboot         = varargin{5};
        clear varargin
        
        % ---------------------------------------------------------------
        % here we work differently because the number of betas etc .. is
        % known only after running the rep_anova
        % Data = [electrode, frames, subjects, conditions]
        % ----------------------------------------------------------------
        % compute the ANOVA
        
        start = 1;
        for electrode = 1:size(data,1)
            
            tmp = squeeze(data(electrode,:,:,:));
            for i=1:size(tmp,3)
                Y(:,:,i) = tmp(:,find(~isnan(tmp(1,:,i))));
                gp = gp_vector(find(~isnan(tmp(1,:,i))));
            end
            
            if ~isempty(gp)
                for i=1:max(gp)
                    ss(i)=sum(gp==i);
                end
            end
            
            if length(unique(ss)) == 1
                fprintf('analyse electrode %g', electrode); disp(' ');
                if start == 1  % use start and not electrode in case electrode 1 would be wrong
                    results = limo_rep_anova(Y,gp,factor_levels,1);
                    if ~isempty(results)
                        % create files to be updated now - design matrix accepted by the user
                        index1 = 0; for i = 1:size(results.Betas,2); index1 = index1 + size(results.Betas{i},1); end
                        Betas = single(NaN(size(data,1), size(data,2), index1));
                        if size(results.F,2) == size(data,2)
                            Repeated_measures_analysis = single(NaN(size(data,1), size(data,2), 2)); % F and p values to store
                        else
                            Repeated_measures_analysis = single(NaN(size(data,1), size(data,2), size(results.F,2), 2)); % F and p values to store
                        end
                        name = sprintf('Repeated_measures_ANOVA');
                        start = 0;
                    else
                        return
                    end
                else
                    results = limo_rep_anova(Y,gp,factor_levels,0);
                end
                clear Y tmp
                
                % update files
                index = 1;
                for i=1:size(results.Betas,2)
                    Betas(electrode,:,index:(index-1)+size(results.Betas{i},1)) = (results.Betas{i})';
                    index = index+size(results.Betas{i},1);
                end
                
                if size(results.F,2) == size(data,2)
                    Repeated_measures_analysis(electrode,:,1) = results.F;
                    Repeated_measures_analysis(electrode,:,2) = results.p;
                else
                    Repeated_measures_analysis(electrode,:,:,1) = results.F;
                    Repeated_measures_analysis(electrode,:,:,2) = results.p;
                end
                
            else
                if start == 0
                    fprintf('electrode %g skipped - groups of non equivalent size due to bad channel', electrode); disp(' ');
                    Betas(electrode,:,:,:) = single(NaN);
                    Repeated_measures_analysis(electrode,:,:,:) = single(NaN);
                    clear Y tmp
                end
            end
        end
        
        [boot_electrode,b]=find(~isnan(Repeated_measures_analysis(:,1,1,1)));
        load LIMO; LIMO.design.X = results.Design; LIMO.Design.name = results.name; save LIMO LIMO;
        save Betas Betas; save ([name],'Repeated_measures_analysis')
        clear LIMO results ss
        
        
        % ----------------------------------------------------------------
        if nboot > 0
            
            % create files to store bootstrap under H1 and H0
            disp('making bootstrap files ...')
            boot_H0_Betas = single(NaN(size(Betas,1), size(Betas,2), size(Betas,3),nboot));
            if numel(size(Repeated_measures_analysis)) == 3
                boot_H0_Repeated_measures_analysis = single(NaN(size(Repeated_measures_analysis,1), size(Repeated_measures_analysis,2), 2, nboot));
            else
                boot_H0_Repeated_measures_analysis = single(NaN(size(Repeated_measures_analysis,1), size(Repeated_measures_analysis,2), size(Repeated_measures_analysis,3), 2, nboot));
            end
            name = sprintf('boot_H0_Repeated_measures_ANOVA');
            clear Repeated_measures_analysis Betas
            
            % the data have to be centered (H0) for each cell
            for i=1:max(gp)
                ss(i)=sum(gp_vector==i);
            end
            
            centered_data = single(NaN(size(data,1),size(data,2),size(data,3),size(data,4)));
            
            gp_index1 = 1;
            nb_conditions = prod(factor_levels);
            for gp=1:length(ss)
                gp_index2 = gp_index1-1+ss(gp);
                for condition=1:nb_conditions
                    centered_data(:,:,gp_index1:gp_index2,condition) = data(:,:,gp_index1:gp_index2,condition);
                end
                gp_index1 = gp_index1+ss(gp);
            end
            
            % create an index to use across all electrodes and frames
            % (different per gp but identical across conditions)
            disp('making random table...')
            index1=1;index2=ss(1);
            boot_index=single(NaN(size(centered_data,3),nboot));
            for cel=1:length(ss)
                if size(centered_data,1)==1
                    chdata(1,:)=squeeze(centered_data(1,1,index1:index2,1)); % used to check data for NaNs - only take 1 condition since NaNs are on electrodes
                else
                    chdata=squeeze(centered_data(:,1,index1:index2,1));
                end
                
                B=1;
                while B~=nboot+1
                    tmp = (ceil(rand(1,ss(cel)).*ss(cel)+index1-1));
                    if length(unique(tmp)) ~= 1 && min(sum(~isnan(chdata(:,tmp-index1+1)),2)) > 1;
                        boot_index(index1:index2,B) = tmp; % works because subjects are stacked by group on dim(3) of data
                        B=B+1;
                    end
                end
                
                if cel<length(ss)
                    index1 = index1+ss(cel);
                    index2 = index2+ss(cel+1);
                end
                clear chdata
            end
            
            % compute under H0
            for B=1:nboot
                
                start = 1;
                for e=1:length(boot_electrode);
                    electrode = boot_electrode(e);
                    tmp2 = squeeze(centered_data(electrode,:,boot_index(:,B),:));
                    for i=1:size(tmp2,3)
                        Y2(:,:,i) = tmp2(:,find(~isnan(tmp2(1,:,i))));
                        gp = gp_vector(find(~isnan(tmp2(1,:,i))));
                    end
                    
                    if ~isempty(gp)
                        for i=1:max(gp)
                            ss(i)=sum(gp==i);
                        end
                    end
                    
                    if length(unique(ss)) == 1
                        fprintf('boot %g analyse electrode %g', B,electrode); disp(' ');
                        results2 = limo_rep_anova(Y2,gp,factor_levels,0);
                        clear tmp2 Y2 % tmp1 Y1
                        
                        % update files
                        index = 1;
                        for i=1:size(results2.Betas,2)
                            boot_H0_Betas(electrode,:,index:(index-1)+size(results2.Betas{i},1),B) = (results2.Betas{i})'; % H0
                            index = index+size(results2.Betas{i},1);
                        end
                        
                        if size(results2.F,2) == size(centered_data,2)
                            boot_H0_Repeated_measures_analysis(electrode,:,1,B) = results2.F; % H0
                            boot_H0_Repeated_measures_analysis(electrode,:,2,B) = results2.p;
                        else
                            boot_H0_Repeated_measures_analysis(electrode,:,:,1,B) = results2.F; % H0
                            boot_H0_Repeated_measures_analysis(electrode,:,:,2,B) = results2.p;
                        end
                    else
                        fprintf('electrode %g skipped - groups of non equivalent size due to bad channel', electrode); disp(' ');
                        clear tmp2 Y2 % tmp1  Y1
                    end
                end
            end
            
            % save boot_H1_Betas boot_H1_Betas; clear boot_H1_Betas
            save ([name],'boot_H0_Repeated_measures_analysis'); clear boot_H0_Repeated_measures_analysis
            save boot_index boot_index % since Yr is too big we can recreate it on the fly for contrasts
        end
        
end



