function limo_1stLevel(designs,opts)

if not(exist('opts','var'))
    opts = struct;
end
def.rootdir = cd;
def.re = '.*.set';
def.rmchan = {'EXG.*'};
def.timewindow = [-Inf Inf];
def.eval = '';
opts = setdef(opts,def);


%%%% first add a few things you might need to the matlab path
global LIMO EEG
cd(opts.rootdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This pattern matches all the files of your subjects.
% see the help of datagrabber.
GRAB = datagrabber(opts.re, 'loadmode','info');
%%% now we have this GRAB structure with one element per subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for i_des = 1:numel(designs)
    cd(opts.rootdir)
    CatStr = [];ContStr = [];
    if not(isempty(designs(i_des).cat))
        for i_cat = 1:numel(designs(i_des).cat)
            CatStr{i_cat} = designs(i_des).cat(i_cat).eval;
        end
    end
    if not(isempty(designs(i_des).cont))
        orthogCont = designs(i_des).cont(1).orthog;
        for i_cont = 1:numel(designs(i_des).cont)
            ContStr{i_cont} = designs(i_des).cont(i_cont).eval;
        end
    end
    
    for isubj = 1:size(GRAB,1)
        LIMO = [];
        EEG = [];
        
        disp(GRAB(isubj).name);
        
        EEG = pop_loadset(GRAB(isubj).name);
        
        if isfield(opts,'reref')
            EEG = pop_reref(EEG,chnb(opts.reref),'keepref','on');
        end
        
        EEG.event = GRAB(isubj).EEG.event;
        % we go to this subject's directory
        sujdir = fullfile(opts.rootdir,GRAB(isubj).suj);
        cd(sujdir);
        
        % Create design directory if necessary
        if not(isdir(designs(i_des).name))
            mkdir(designs(i_des).name);
        end
        cd(designs(i_des).name)
        disp(designs(i_des).name)
        delete Betas.* Betas_* con_* ess_* Condition_effect.* Continuous.* R2.* Res.* Yhat.*mat Yr.*mat
        %%% remove unwanted channels.
        EEG = pop_select( EEG,'nochannel',chnb(opts.rmchan));
        
        % select events of interest
        goodevents = [];
        for i = 1:numel(designs(i_des).events_of_interest)
            goodevents = [goodevents find(cellstr2num({EEG.event.type}) == designs(i_des).events_of_interest(i))];
        end
        goodevents = sort(goodevents);
        event = EEG.event(goodevents);
        for i_ev = 1:numel(event)
            event(i_ev).type = num2str(event(i_ev).type);
        end
        
        if isfield(designs,'eval') && not(isempty(designs.eval))
            eval(designs(i_des).eval);
        end
        
        %%%% putting the Categorical conditions together
        LIMOCat = [];
        if not(isempty(designs(i_des).cat))
            LIMOCat = zeros(1,numel(event));
            for i = 1:numel(CatStr)
                LIMOCat(eval(CatStr{i})) = i;
            end
        end
        LIMOCat(LIMOCat == 0) = max(LIMOCat)+1;
        if max(LIMOCat(:)) > numel(CatStr)
            disp('Some trials are unmodelled. An additional regressor has been added in last position.')
        end
        save LIMOCat LIMOCat
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% putting the Continuous conditions together
        LIMOCont = [];
        if not(isempty(designs(i_des).cont))
            for i_cont = 1:numel(ContStr)
                tmp = eval(ContStr{i_cont});
                if size(tmp,2) > size(tmp,1)
                    tmp = tmp';
                end
                LIMOCont = [LIMOCont tmp];
            end
            LIMOCont = zscore(LIMOCont); % zscoreing
            % add interaction?
            if numel(designs(i_des).cont(1).interact) == 1 && designs(i_des).cont(1).interact
                if numel(designs(i_des).cont) ~= 2
                    error('check interactions')
                end
                LIMOCont(:,end+1) = LIMOCont(:,1).*LIMOCont(:,2);
            elseif numel(designs(i_des).cont(1).interact) > 1
                [idx jdx] = find(designs(i_des).cont(1).interact);
                for i = 1:numel(idx)
                    LIMOCont(:,end+1) = LIMOCont(:,idx(i)).*LIMOCont(:,jdx(i));
                end
            end
            if size(LIMOCont,2) > 1
                disp('Correlation between regressors')
                [rho pval] = corr(LIMOCont)
                % orthogonalizing them if required
                if orthogCont 
                    LIMOCont = limo_orth(LIMOCont);
                end
            end
        end
        save LIMOCont LIMOCont
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [timesofinterest tvals] = timepts(opts.timewindow,EEG.times);
        
        % initialize LIMO structure
        LIMO.data.data_dir          = sujdir;
        LIMO.dir                    = cd;
        LIMO.data.data              = GRAB(isubj).name;
        LIMO.data.chanlocs          = EEG.chanlocs;
        LIMO.data.start             = tvals(1)/1000;
        LIMO.data.end               = tvals(end)/1000;
        LIMO.data.sampling_rate     = EEG.srate;
        LIMO.data.Cat               = LIMOCat;
        LIMO.data.Cont              = LIMOCont;
        LIMO.data.trim1             = timesofinterest(1);
        LIMO.data.trim2             = timesofinterest(end);
        LIMO.Method                 = 1;
        LIMO.design.zscore          = 0;% I did it above.
        if isfield('intercept',designs)
            LIMO.design.intercept       = designs(i_des).intercept;
        else
            LIMO.design.intercept       = 1;
        end
        LIMO.Bootstrap              = 0;
        LIMO.Level                  = 1;
        
        
        % create design matrix
        Y = EEG.data(:,LIMO.data.trim1:LIMO.data.trim2,:);
        
        [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_continuous, ...
            LIMO.design.nb_items, LIMO.design.conditions] = ...
            limo_design_matrix(Y, LIMO);
        %%% limo_design_matrix does all sorts of crasy things like sorting the
        %%% data for some unknown reason.
        
        % display it for the first subject, just to check.
        if isubj == 1
            figure('Name',['LIMO design matrix ' EEG.subject],'Color','w','NumberTitle','off')
            Xdisplay = LIMO.design.X;
            if   LIMO.design.nb_continuous
                REGdisplay = LIMO.design.X(:, LIMO.design.nb_conditions+1:size(LIMO.design.X,2)-LIMO.design.intercept);
                REGdisplay = REGdisplay + max(abs(min(REGdisplay)));
                Xdisplay(:, LIMO.design.nb_conditions+1:size(LIMO.design.X,2)-LIMO.design.intercept) = REGdisplay ./ max(max(REGdisplay));
            end
            imagesc(Xdisplay); colormap('gray'); drawnow;
            title('Design matrix'); xlabel('regressors');ylabel('trials');
            set(gca,'XTick',1:size( LIMO.design.X,2))
            title(GRAB(isubj).suj,'interpreter','none');
            drawnow
        end
        
        % defining the type of analysis that will be run, based on the
        % regressors that have been defined.
        if LIMO.design.nb_conditions > 0 && LIMO.design.nb_continuous == 0
            if LIMO.design.nb_conditions == 2
                LIMO.design.name  = sprintf('Categorical: T-test i.e. %g conditions',LIMO.design.nb_conditions);
            else
                LIMO.design.name  = sprintf('Categorical: ANOVA with %g conditions',LIMO.design.nb_conditions);
            end
        elseif LIMO.design.nb_conditions == 0 && LIMO.design.nb_continuous > 0
            if LIMO.design.nb_continuous == 1
                LIMO.design.name  = sprintf('Continuous: Simple Regression');
            else
                LIMO.design.name  = sprintf('Continuous: Multiple Regression with %g continuous variables',LIMO.design.nb_continuous);
            end
            
        elseif LIMO.design.nb_conditions > 0 && LIMO.design.nb_continuous > 0
            LIMO.design.name      = sprintf('AnCOVA with %g conditions and %g continuous variable(s)',LIMO.design.nb_conditions,LIMO.design.nb_continuous);
        else
            error('Don''t know what to do.');
        end
        
        disp('design matrix done ...')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % now running the analysis
        % --------- load files created by limo_design_matrix ------------------
        load Yr; load Yhat; load Res; load R2; load Betas;
        
        if LIMO.design.nb_conditions ~=0
            load Condition_effect
        end
        
        if LIMO.design.nb_continuous ~=0
            load Continuous
        end
        
        % -------------- loop the analysis electrode per electrode
        % ---------------  GLM Part ------------------------------
        for electrode = 1:size(Yr,1)
            warning off;
            if electrode>1
                for i = 1:22
                    fprintf('\b');
                end
            end
            fprintf('analyzing electrode %02g',electrode);
            if electrode == size(Yr,1)
                fprintf('\n');
            end
            %%%%% this is the one who's doing the analysis.
            model = limo_glm(squeeze(Yr(electrode,:,:))',LIMO); warning on;
            
            
            % update the LIMO degrees of freedom (do it only once)
            if electrode == 1
                LIMO.model.model_df = model.df;
                if LIMO.design.nb_conditions ~=0 % everything is in model
                    LIMO.model.conditions_df  = model.univariate.conditions.df;
                end
                if LIMO.design.nb_continuous ~=0
                    LIMO.model.continuous_df  = model.univariate.continuous.df;
                end
            end
            
            % update the files to be stored on the disk
            fitted_data = LIMO.design.X*model.betas;
            % this is all the data that has been accounted for by the model.
            Yhat(electrode,:,:) = fitted_data';
            % Below are the residuals (data that has not been accounted for by
            % the model).
            Res(electrode,:,:)  = squeeze(Yr(electrode,:,:)) - fitted_data'; clear fitted_data
            % R2 is the %variance accounted for by the model. We store the F and associated p value
            R2(electrode,:,1) = model.R2_univariate; R2(electrode,:,2) = model.F; R2(electrode,:,3) = model.p;
            % betas are the regression coefficients.
            Betas(electrode,:,:) = model.betas';
            
            if LIMO.design.nb_conditions ~=0
                % this is the variance explained by the model's categorical
                % variables (all of them together)
                Condition_effect(electrode,:,1) = model.univariate.conditions.F;
                Condition_effect(electrode,:,2) = model.univariate.conditions.p;
            end
            
            if LIMO.design.nb_continuous ~=0
                for i=1:LIMO.design.nb_continuous
                    % this is the variance explained by the model's continuous
                    % variables (one by one)
                    Continuous(electrode,:,i,1) = model.univariate.continuous.F(i,:);
                    Continuous(electrode,:,i,2) = model.univariate.continuous.p(i,:);
                end
            end
        end
        
        % save data on the disk
        save Yhat Yhat;
        save Res Res;
        save Betas Betas;
        save R2 R2;
        save LIMO LIMO;
        if LIMO.design.nb_conditions ~=0
            save Condition_effect Condition_effect
        end
        
        if LIMO.design.nb_continuous ~=0
            save Continuous Continuous
        end
        clear file electrode filename model reg dir i
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% saving in EEG format.
        % we are creating EEG sets that contain only one epoch each. This is
        % considered continuous data by eeglab and will prevent using certain
        % plotting tools later on. Handle with care.
        
        % first creating an empty EEG set
        EEG = eeg_emptyset;
        EEG.chanlocs = LIMO.data.chanlocs;
        EEG.srate = LIMO.data.sampling_rate;
        EEG.xmin = LIMO.data.start;
        EEGi = EEG;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if LIMO.design.nb_conditions
            EEG.setname = 'Condition_effect (F values)';
            EEG.comments = sprintf('Condition_effect: %s',fullfile(cd,'Condition_effect.mat'));
            EEG.datfile = fullfile(cd,'Condition_effect.mat');
            
            EEG = mat2eeglab(Condition_effect(:,:,1),EEG);
            pop_saveset(EEG, 'filename', ['Condition_effect.set'],'filepath',cd);
            clear Condition_effect
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if LIMO.design.nb_continuous
            EEG.setname = 'Continuous (F values)';
            EEG.comments = sprintf('Continuous: %s',fullfile(cd,'Continuous.mat'));
            EEG.datfile = fullfile(cd,'Continuous.mat');
            
            EEG = mat2eeglab(Continuous(:,:,1),EEG);
            pop_saveset(EEG, 'filename', ['Continuous.set'],'filepath',cd);
            clear Continuous
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iB = 1:size(Betas,3)
            EEG = EEGi;
            EEG.setname = ['Betas (' num2str(iB) ')'];
            EEG.comments = sprintf('Betas: %s',fullfile(cd,'Betas.mat'));
            
            EEG = mat2eeglab(Betas(:,:,iB),EEGi);
            pop_saveset(EEG, 'filename', ['Betas_' num2str(iB,'%02d') '.set'],'filepath',cd);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % here we are also catenating all the Betas as epochs in one dataset.
        % could be useful to look at ERP images.
        EEG = EEGi;
        EEG.setname = ['Betas (all)'];
        EEG.comments = sprintf('Betas: %s',fullfile(cd,'Betas.mat'));
        
        EEG = mat2eeglab(Betas,EEGi);
        pop_saveset(EEG, 'filename', ['Betas.set'],'filepath',cd);
        clear Betas
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EEG = EEGi;
        EEG.setname = ['Residuals'];
        EEG.comments = sprintf('Residuals: %s',fullfile(cd,'Res.mat'));
        
        EEG = mat2eeglab(Res,EEGi);
        pop_saveset(EEG, 'filename', ['Res.set'],'filepath',cd);
        clear Res
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EEG = EEGi;
        EEG.setname = ['Rsquare'];
        EEG.comments = sprintf('Rsquare: %s',fullfile(cd,'R2.mat'));
        
        EEG = mat2eeglab(R2(:,:,1),EEGi);
        pop_saveset(EEG, 'filename', ['R2.set'],'filepath',cd);
        clear R2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
disp done.
cd(opts.rootdir)
clear GRAB EEG

function s = setdef(s,d)
% s = setdef(s,d)
% Merges the two structures s and d recursively.
% Adding the default field values from d into s when not present or empty.

if isstruct(s) && not(isempty(s))
    fields = fieldnames(d);
    for i_f = 1:numel(fields)
        if isfield(s,fields{i_f})
            s.(fields{i_f}) = setdef(s.(fields{i_f}),d.(fields{i_f}));
        else
            s.(fields{i_f}) = d.(fields{i_f});
        end
    end
elseif not(isempty(s))
    s = s;
elseif isempty(s);
    s = d;
end

