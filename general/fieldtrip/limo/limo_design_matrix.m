function [X, nb_conditions, nb_continuous, nb_items, conditions] = limo_design_matrix(varargin)

% LIMO_DESIGN_MATRIX - called once by LIMO_EEG.m to create the design matrix X that is
% stored in the LIMO.mat file and called for all analyses.
%
% FORMAT:
% [X,nb_conditions,nb_continuous,nb_items] = limo_design_matrix(Y,LIMO)
% [X,nb_conditions,nb_continuous,nb_items] = limo_design_matrix(Y,Cat,Cont,directory,zscoring)
%
% INPUTS: 
%   Y             = EEG data with format electrodes x frames x trials/subjects
%   Cat           = vector describing the different conditions (if no conditions Cat = 0)
%   Cont          = matrix describing the different covariates (if no covariates Cont = 0)
%   directory     = path of folder where the outputs will be saved (see below)
%   zscoring      = [0/1] - if 1 continuous regressors are zscored, 
%                          which means that the betas coefficients will have units
%                          micro-volts per std of the predictor variable.
%   intercept     = [0/1] include an intercept or not in the design matrix
%   LIMO          = = structure that contains the above information (except Y)
%
% OUTPUTS: 
%   X             = 2 dimensional matrix that describes the experiments' events
%   nb_conditions = returns the number of conditions / groups
%   nb_continuous = returns the number of covariates / regressors
%   nb_items      = returns the number of trials per condition (allows unbalanced designs)
%   conditions    = conditions present in the data
%
%   These outputs are written to disk in DIRECTORY and populated latter:
%
%   Yr.mat = the EEG data reorganized, grouped by conditions (if Cat ~=0)
%   Yhat.mat = the predicted data (same size as Y)
%   Res.mat = the residual (non modeled) data (same size as Y)
%   R2.mat = the model fit (size=[1,frame,electrode])
%   Beta_00X.mat = the beta values (size=[1,frame,electrode])
%
% See also LIMO_EEG LIMO_GLM
%
% Cyril Pernet / Guillaume Rousselet v4 27/04/2009 
% -----------------------------
%  Copyright (C) LIMO Team 2010


%% varagin stuffs

if length(varargin)==2
    Y         = varargin{1}; 
    Cat       = varargin{2}.data.Cat;
    Cont      = varargin{2}.data.Cont;
    directory = varargin{2}.dir;
    zscoring  = varargin{2}.design.zscore;
    intercept = varargin{2}.design.intercept;
elseif length(varargin)==5
    Y         = varargin{1}; 
    Cat       = varargin{2};
    Cont      = varargin{3};
    directory = varargin{4};
    zscoring  = varargin{5};
else
    error('varargin error')
end


%% Basic checking

test = size(Cat) + size(Cont);
if test ==[2 2]
    errordlg('no regressors selected','file selection error'); return
end

% check Cat
[l,w]=size(Cat);
if l == 1 && w == size(Y,3)
    Cat = Cat';
    disp('Cat has been transposed')
end

% check Cont
[l,w]=size(Cont);
if w == size(Y,3)
    Cont = Cont';
    disp('Cont has been transposed')
end

% overall dimensions check
if isempty(Cat)
    if size(Y,3) ~= length(Cont)
        error('The number of trials and the covariate(s) length are different')
    end
elseif isempty(Cont)
    if size(Y,3) ~= length(Cat)
        error('The number of trials and the number of events are different size')
    end
else % cat and cont ~= 0
    if size(Y,3) ~= length(Cont)
        error('The number of trials and the covariate(s) length are different')
    elseif size(Y,3) ~= length(Cat)
        error('The number of trials and the number of events are different size')
    elseif length(Cat) ~= length(Cont)
        error('The number of events and the covariate(s) length are different')
    end
end

% additional checking for regressions
if Cat == 0

    if size(Cont,2)+1 >= size(Y,3) 
        error('there are too many regressors for the number of trials, reduce the model size')
        % one cannot compute any statistics if size(Y,3) > size(Cont,2)+1
        % because of the dof - could do up to some levels using a Tikhonov
        % regularization but not tested for now
    end

    if size(Cont,2) >1
        if det(Cont'*Cont) == 0
            errordlg('the regression matrix is singular, 1 or more predictors are function to each other, please consider another model', ...
                'Matrix definition error'); return
        elseif cond(Cont'*Cont,1) == Inf
            errordlg('the regression matrix is close to singular, please consider another model','Matrix definition error'); return
        end

    end
end

% cd(directory);


%% Make the design matrix and create files


conditions = unique(Cat);
nb_conditions = numel(conditions); % stands for the categorical variable
nb_continuous = size(Cont,2); % stands for the continuous regressors
nb_items      = 0; % stands for sample size per condition
if nb_continuous
    [l,w]=size(Cont);
    if l~=size(Y,3)
        errordlg('The covariate(s) must be the same length as your dependant variable(s)')
        error('Please retry changing the 3rd argument');
    end
    Continuous = zeros(size(Y,1),size(Y,2),nb_continuous,2,'single');
    save Continuous Continuous; clear Continuous
end
if nb_conditions
    [l,w]=size(Cat);
    if w ~=1
        errordlg('The categorical regressor must be a vector, i.e. 1 column')
        error('please retry changing the 2nd argument');
    elseif l~=size(Y,3)
        errordlg('The categorical regressor must be the same length as your dependant variables')
        error('please retry changing the 2nd argument');
    end
    Condition_effect = zeros(size(Y,1),size(Y,2),2,'single'); % dim 3 = F/p values
    save Condition_effect Condition_effect; clear Condition_effect;
end

disp('Creating the design matrix and data files ...')
X = zeros(size(Y,3),nb_conditions+nb_continuous+intercept);
for i = 1:nb_conditions
    X(:,i) = double(Cat == conditions(i));
end
if zscoring == 1
    Cont = zscore(Cont);
end
j = 1;
for i = nb_conditions+1:nb_conditions+nb_continuous
    X(:,i) = Cont(:,j);
    j = j+1;
end
if intercept
    X(:,end) = 1;
end
nb_items = sum(X(:,1:nb_conditions));
% if not(isempty(Cat))
%     [dum order] = sort(Cat);
% else
%     order = 1:size(X,1);
% end
% X = X(order,:);
Yr = Y;save Yr Yr; clear Yr%(:,:,order);                                      save Yr Yr; clear Yr
Yhat  = zeros(size(Y),'single');                        save Yhat Yhat; clear Yhat
Res   = zeros(size(Y),'single');                        save Res Res; clear Res
Betas = zeros(size(Y,1),size(Y,2),size(X,2),'single');  save Betas Betas; clear Betas
R2    = zeros(size(Y,1),size(Y,2),3,'single');          save R2 R2; clear R2

% % % figure
% figure('Name','LIMO design','Color','w','NumberTitle','off')
% Xdisplay = X; 
% if  nb_continuous
%     REGdisplay = X(:,nb_conditions+1:size(X,2)-1); 
%     REGdisplay = REGdisplay + max(abs(min(REGdisplay)));
%     Xdisplay(:,nb_conditions+1:size(X,2)-1) = REGdisplay ./ max(max(REGdisplay));
% end
% imagesc(Xdisplay); colormap('gray'); drawnow; 
% title('Design matrix'); xlabel('regressors');ylabel('trials');
% set(gca,'XTick',1:size(X,2))
% % 

