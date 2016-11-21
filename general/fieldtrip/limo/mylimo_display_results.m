function mylimo_display_results(Type,FileName,PathName,p,MCC,LIMO,varargin)

% This function displays various results
% The arguments specify cases for the
% different kind of figures, thresholds etc ..
%
% FORMAT:
% limo_display_results(Type,FileName,PathName,p,LIMO)
%
% INPUTS:
%   Type      = type of images/plot to do
%   Filename  = Name of the file to image
%   PathName  = Path of the file to image
%   p         = threshold p value from the limo_results GUI, e.g. 0.05
%   MCC       = Multiple Comparison technique
%   LIMO      = LIMO structure
%
% Types of images
%  1 - 2D images with a intensity plotted as function of time (x) and electrodes (y)
%  2 - topographic plot a la eeglab
%  3 - plot the ERP data (original or modeled)
%
% Cyril Pernet, Guillaume Rousselet v3 06-05-2009
% Carl Gaspar 03-09-2009 - fixed some axis issues for 3D plots (see subfunction time_vect_)
% Cyril P. v4 09-09-2009 allows random effect results to be displayed (+ some clean up)
% Cyril P. v5. 10-03-2010 split the whole file into 2 parts based on LIMO.level (1 or 2)
% Guillaume Rousselet v4 06-07-2010 added the max(T)/max(F) and cluster stats for random effect
% Cyril Pernet v4 16-05-2010 fixed the random effect to automatically load bootstrap and get the neighbouring matrix for clusters
% -----------------------------
%  Copyright (C) LIMO Team 2010
 
 
% cd(PathName)
[p f e] = fileparts(FileName);
switch e
    case '.mat'
        load (FileName);
    case '.set'
        EEG = pop_loadset(FileName,PathName);
end
% -------------------------------------------------------------------------
% -------------------      LEVEL 1     ------------------------------------
% -------------------  SINGLE SUBJECT  ------------------------------------
% -------------------------------------------------------------------------
if LIMO.Level == 1
    
    switch Type
        
        case{1}
            
            %--------------------------
            % imagesc of the results
            %--------------------------
            
            % univariate results from 1st level analysis
            % ------------------------------------------
            
            if strcmp(FileName,'R2.mat')
                mask = squeeze(R2(:,:,3)) < p;
                if sum(sum(sum(mask))) == 0
                    warndlg('no values under threshold','no significant effect');
                else
                    toplot=squeeze(R2(:,:,1));
                    scale = toplot.*mask; scale(scale==0)=NaN;
                    mytitle = 'R^2 Values';
                    assignin('base','Plotted_data',squeeze(R2(:,:,1)).*mask)
                end
                
            elseif strcmp(FileName(1:end-6),'Condition_effect')
                mask = squeeze(Condition_effect(:,:,2)) < p;
                if sum(sum(sum(mask))) == 0
                    warndlg('no values under threshold','no significant effect');
                else
                    toplot = squeeze(Condition_effect(:,:,1));
                    scale = toplot.*mask; scale(scale==0)=NaN;
                    mytitle = sprintf('Condition effect %s - F Values',FileName(end-4));
                    tmp_name = sprintf('Plotted_data_condition_%s',FileName(end-4)');
                    assignin('base',[tmp_name],squeeze(Condition_effect(:,:,1)).*mask)
                end
                
            elseif strcmp(FileName(1:end-6),'Covariate_effect') 
                mask = squeeze(Covariate_effect(:,:,2)) < p;
                if sum(sum(sum(mask))) == 0
                    warndlg('no values under threshold','no significant effect');
                else
                    toplot = squeeze(Covariate_effect(:,:,2));
                    scale = toplot.*mask; scale(scale==0)=NaN;
                    mytitle = sprintf('Covariate %s - F values',FileName(end-4));
                    tmp_name = sprintf('Plotted_data_continuous_variable_%s',FileName(end-4)');
                    assignin('base',[tmp_name],squeeze(Covariate_effect(:,:,2)).*mask)
                end
                
            elseif strcmp(FileName,'Partial_coef.mat')
                for i=1:size(Partial_coef,3)
                    mask = squeeze(Partial_coef(:,:,i,3)) < p;
                    if sum(mask(:)) == 0
                        warndlg('no values under threshold','no significant effect');
                    else
                        toplot = squeeze(Partial_coef(:,:,i,1));
                        scale = toplot.*mask; scale(scale==0)=NaN;
                        mytitle = sprintf('Partial coef R2 values variable %g',i);
                        tmp_name = sprintf('Plotted_data_partial_coef_%g',i);
                        assignin('base',[tmp_name],squeeze(Partial_coef(:,:,i,1)).*mask)
                    end
                end
                
            elseif strcmp(FileName(1:4),'con_')
                mask = squeeze(con(:,:,3)) < p;
                if sum(sum(sum(mask))) == 0
                    warndlg('no values under threshold','no significant effect');
                else
                    toplot = squeeze(con(:,:,2));
                    scale = toplot.*mask; scale(scale==0)=NaN;
                    mytitle = ['Contrast ',[FileName(5:end-4)],' -- T values'];
                    assignin('base','Plotted_data',squeeze(con(:,:,2)).*mask)
                end
                
            elseif strcmp(FileName(1:4),'ess_')
                mask = squeeze(ess(:,:,3)) < p;
                if sum(sum(sum(mask))) == 0
                    warndlg('no values under threshold','no significant effect');
                else
                    toplot = squeeze(ess(:,:,end-1));
                    scale = toplot.*mask; scale(scale==0)=NaN;
                    mytitle = ['Contrast ',[FileName(5:end-4)],' -- F values'];
                    assignin('base','Plotted_data',squeeze(ess(:,:,end-1)).*mask)
                end
                
            else
                disp('file not supported');
                return
            end
            
            figure;set(gcf,'Color','w');
            ax(1) = subplot(3,3,[1 2 4 5 7 8]);
            timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(toplot,2));
            ratio = (LIMO.data.end*1000 - LIMO.data.start*1000) / size(toplot,2);
            if LIMO.data.start < 0
                frame_zeros = round(abs(LIMO.data.start*1000) / ratio);
            end
            imagesc(timevect,1:size(toplot,1),toplot.*mask);
            v = max(toplot(:)); [e,f]=find(scale==v);
            colormap(hot);
            try
                caxis([min(min(scale)), max(max(scale))]);
            end
            color_images_(scale,LIMO);
            title(mytitle,'Fontsize',18)
            
            ax(2) = subplot(3,3,6);
            topoplot(toplot(:,f),LIMO.data.chanlocs);
            title(['topoplot @' num2str(round(timevect(f))) 'ms'],'FontSize',14)
            
            ax(3) = subplot(3,3,9);
            plot(toplot(e,:),'LineWidth',3); grid on
            mytitle2 = sprintf('time course @ \n electrode %s', LIMO.data.chanlocs(e).labels);
            title(mytitle2,'FontSize',14)
            
            update = 0;
            while update ==0
                try
                    [x,y,button]=ginput(1);
                catch
                    update =1; break
                end
                if button > 1
                    update = 1;
                end
                clickedAx = gca;
                if clickedAx ~=ax(1)
                    disp('right click to exit')
                else
                    frame = round(x / ratio);
                    if frame < 0
                        frame = frame_zeros + frame;
                    end
                    subplot(3,3,6,'replace');
                    topoplot(toplot(:,frame),LIMO.data.chanlocs);
                    title(['topoplot @ ' num2str(round(x)) 'ms'],'FontSize',14)
                    
                    y = round(y);
                    subplot(3,3,9,'replace'); 
                    plot(toplot(y,:),'LineWidth',3); grid on
                    try
                        mytitle2 = sprintf('time course @ \n electrode %s', LIMO.data.expected_chanlocs(y).labels);
                    catch ME
                        mytitle2 = sprintf('time course @ \n electrode %s', LIMO.data.chanlocs(y).labels);
                    end
                    title(mytitle2,'FontSize',14)
                end
            end
            
            
        case{2}
            
            %--------------------------
            % topoplot
            %--------------------------
            
            % univariate results from 1st level analysis
            % ------------------------------------------
            
            timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
            
            if strcmp(FileName,'R2.mat')
                EEG.data = squeeze(R2(:,:,1));
                EEG.setname = 'R2 values';
                EEG.pnts = size(EEG.data,2);
                EEG.xmin = LIMO.data.start;
                EEG.xmax = LIMO.data.end;
                EEG.times = timevect;
                EEG.trials = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                EEG.nbchan = size(EEG.data,1);
                pop_topoplot(EEG);
                assignin('base','Plotted_data',squeeze(R2(:,:,1)))
                
            elseif strcmp(FileName(1:end-6),'Condition_effect')
                EEG.data = squeeze(Condition_effect(:,:,1));
                EEG.setname = sprintf('Condition %s - F values',FileName(end-4));
                EEG.pnts = size(EEG.data,2);
                EEG.xmin = LIMO.data.start;
                EEG.xmax = LIMO.data.end;
                EEG.times = timevect;
                EEG.trials = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                EEG.nbchan = size(EEG.data,1);
                pop_topoplot(EEG);
                tmp_name = sprintf('Plotted_data_condition_%s',FileName(end-4));
                assignin('base',[tmp_name],squeeze(Condition_effect(:,:,1)))
                
            elseif strcmp(FileName(1:end-6),'Covariate_effect')
                EEG.data = squeeze(Covariate_effect(:,:,1));
                EEG.setname = sprintf('Covariate %s - F values',FileName(end-4));
                EEG.pnts = size(EEG.data,2);
                EEG.xmin = LIMO.data.start;
                EEG.xmax = LIMO.data.end;
                EEG.times = timevect;
                EEG.trials = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                EEG.nbchan = size(EEG.data,1);
                pop_topoplot(EEG);
                tmp_name = sprintf('Plotted_data_covariate_%s',FileName(end-4));
                assignin('base',[tmp_name],squeeze(Continuous(:,:,b,1)))
                
            elseif strcmp(FileName,'Partial_coef.mat')
                regressor = str2num(cell2mat(inputdlg('which regressor to plot (e.g. 1:3)','Plotting option')));
                if max(regressor) > size(Partial_coef,3); errordlg('error in regressor number'); return; end
                for b = regressor
                    EEG.data = squeeze(Partial_coef(:,:,b,1));
                    EEG.setname = sprintf('Partial coef R2 values variable %g',b);
                    EEG.pnts = size(EEG.data,2);
                    EEG.xmin = LIMO.data.start;
                    EEG.xmax = LIMO.data.end;
                    EEG.times = timevect;
                    EEG.trials = 1;
                    EEG.chanlocs = LIMO.data.chanlocs;
                    EEG.nbchan = size(EEG.data,1);
                    pop_topoplot(EEG);
                    tmp_name = sprintf('Plotted_data_partial_coef_%g',i);
                    assignin('base',[tmp_name],squeeze(Partial_coef(:,:,b,1)))
                end
                
            elseif strcmp(FileName(1:4),'con_')
                if ~exist('EEG','var')
                    toplot = squeeze(con(:,:,2));
                    EEG.setname = ['Contrast ',[FileName(5:end-4)],' -- T values'];
                    EEG.data = toplot;
                    EEG.pnts = size(EEG.data,2);
                    EEG.xmin = LIMO.data.start;
                    EEG.xmax = LIMO.data.end;
                    EEG.times = timevect;
                    EEG.trials = 1;
                    EEG.chanlocs = LIMO.data.chanlocs;
                    EEG.nbchan = size(EEG.data,1);
                end
                pop_topoplot(EEG);
                assignin('base','Plotted_data',squeeze(con(:,:,2)))
                
            elseif strcmp(FileName(1:4),'ess_')
                toplot = squeeze(ess(:,:,end-1));
                EEG.setname = ['Contrast ',[FileName(5:end-4)],' -- F values'];
                EEG.data = toplot;
                EEG.pnts = size(EEG.data,2);
                EEG.xmin = LIMO.data.start;
                EEG.xmax = LIMO.data.end;
                EEG.times = timevect;
                EEG.trials = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                EEG.nbchan = size(EEG.data,1);
                pop_topoplot(EEG);
                assignin('base','Plotted_data',squeeze(ess(:,:,end-1)))
            else
                disp('file not supported');
            end
            
            
            
        case{3}
            
            %--------------------------
            % ERP
            %--------------------------
            
            % which ERP to make
            regressor = inputdlg('which regressor to plot (e.g. 1:3 or all for factorial)','Plotting option');
            if isempty(regressor)
                return;
            end
            timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
          
            
            % Factorial design - multiply cell to make new X or if simple
            % design do all x in X
            if strcmp(regressor,'all')
                regressor = 1:prod(LIMO.design.nb_conditions);
                electrode = inputdlg('which electrode to plot','Plotting option');
                try
                    electrode = eval(cell2mat(electrode));
                catch ME
                    load R2; [v,e] = max(R2(:,:,1)); [v,c]=max(v); electrode = e(c);
                    [v,electrode] = max(max(mean(R2,3),[],2)); % plot at the max of Yr
                end
                
                if size(electrode) > 1
                    errordlg('invalid electrode number')
                elseif electrode > size(LIMO.data.chanlocs,2)
                    errordlg('invalid electrode number')
                else
                    
                    load Yr;  stop = size(Yr,2);
                    figure;set(gcf,'Color','w')
                    
                    X = LIMO.design.X(:,1:prod(LIMO.design.nb_conditions));
                    if length(LIMO.design.nb_conditions) > 1 % make a new design matrix per conditions
                        % get each factor
                        index = 1;
                        for f = 1:length(LIMO.design.nb_conditions)
                            if f == 1
                                factor{f} = X(:, index:LIMO.design.nb_conditions(f));
                                index = index + LIMO.design.nb_conditions(f);
                            else
                                factor{f} = X(:, index:(LIMO.design.nb_conditions(f)+index-1));
                                index = index + LIMO.design.nb_conditions(f);
                            end
                        end
                        % multiply factors
                        clear tmp; index = 1;
                        combination = nchoosek([1:length(LIMO.design.nb_conditions)],length(LIMO.design.nb_conditions));
                        for f = 1:(max(combination)-1)
                            update = [];
                            for c = 1:size(factor{f},2)
                                col = repmat(factor{f}(:,c),1,size(factor{f+1},2));
                                tmp{index} = col.* factor{f+1};
                                update{c} = index;
                                index = index + 1;
                            end
                            v = cell2mat(update);
                            tmp2 = [];
                            for c=1:length(v)
                                tmp2 = [tmp2 tmp{v(c)}];
                            end
                            factor{f+1} = tmp2;
                        end
                        clear X; X = [factor{f+1} ones(size(col,1),1)];
                    end
                    clear tmp tmp2
                    
                    for i=regressor
                        Reg = squeeze(X(:,i));
                        tmp = squeeze(Yr(electrode,:,find(Reg)));
                        ERP(i,:) = squeeze(mean(tmp,2));   % categorical variable 2D plot
                        clear tmp;
                    end
                end
                
                plot(timevect,ERP','LineWidth',1.5);
                mytitle = sprintf('ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                assignin('base','Plotted_data', ERP')
                grid on; axis tight; title(mytitle,'FontSize',19); drawnow;
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                ylabel('Amplitude in {\mu}V','FontSize',16)
                xlabel('Time in ms','FontSize',16)
                
                
            else % just use the current design X
                try
                    regressor = eval(cell2mat(regressor));
                catch ME
                    return
                end
                
                if max(regressor) > size(LIMO.design.X,2)
                    errordlg('invalid regressor number')
                else
                    extra = questdlg('Plotting ERP','ERP Options','Original','Modeled','Adjusted','Adjusted');
                    if isempty(extra)
                        return;
                    end
                    electrode = inputdlg('which electrode to plot','Plotting option');
                    try
                        electrode = eval(cell2mat(electrode));
                    catch ME
                        load R2; [v,e] = max(R2(:,:,1)); [v,c]=max(v); electrode = e(c);
                        % [v,electrode] = max(max(mean(Yr,3),[],2)); % plot at the max of Yr
                    end
                    
                    if size(electrode) > 1
                        errordlg('invalid electrode number')
                    elseif electrode > size(LIMO.data.chanlocs,2)
                        errordlg('invalid electrode number')
                    else
                        
                        if length(regressor) ~= 1 && (max(regressor) > sum(LIMO.design.nb_conditions))
                            errordlg('you cannot plot categorical and continuous variables at the same time','error');
                            return
                        end
                        load Betas ; load Yr;  stop = size(Betas,2);
                        
                        figure;set(gcf,'Color','w')
                        
                        if strcmp(extra,'Modeled')
                            load Res; subplot(5,1,5);
                            plot(timevect,mean(squeeze(Res(electrode,1:length(timevect),:)),2));
                            grid on; axis tight; title('Model Residuals','FontSize',10); xlabel('Time in ms','FontSize',16)
                        end
                        
                        
                        for i=regressor
                            for frame=1:stop
                                if strcmp(extra,'Modeled')
                                    subplot(5,1,1:4);
                                    % Yhat(frame,:) = LIMO.design.X*squeeze(Betas(electrode,frame,:));% contains all trials
                                    % tmp(frame,:) = squeeze(LIMO.design.X(:,i)).*squeeze(Yhat(frame,:)') ;
                                    Reg = squeeze(LIMO.design.X(:,[i end]));
                                    tmp(frame,:) = Reg(:,1).*squeeze(Betas(electrode,frame,i)) + Reg(:,2).*squeeze(Betas(electrode,frame,end)); % model of param + model of cst
                                else % strcmp(extra,'Original') and strcmp(extra,'Adjusted')
                                    Reg = squeeze(LIMO.design.X(:,i));
                                    tmp(frame,:) = (Reg(find(Reg)).*squeeze(Yr(electrode,frame,find(Reg)))) ;
                                end
                            end
                            
                            if i <= sum(LIMO.design.nb_conditions)
                                ERP(i,:) = squeeze(mean(tmp,2));   % categorical variable 2D plot
                            else
                                [data,index]=sort(LIMO.design.X(:,i));  % continuous variable 3D plot
                                ERP(i,:,:) = tmp(:,index);
                            end
                            clear tmp;
                            
                            if strcmp(extra,'Adjusted')
                                reduced = eye(size(LIMO.design.X,2)-1);
                                of_no_interest = setdiff([1:size(reduced,2)],regressor);
                                for frame=1:stop
                                    tmp2(frame,:) = mean(squeeze(LIMO.design.X(:,of_no_interest)).*repmat(squeeze(Yr(electrode,frame,:)),1,length(of_no_interest)),2);
                                end
                                ERP2(i,:) = squeeze(mean(tmp2,2)); subplot(5,1,5);
                                plot(timevect,ERP2','--','LineWidth',1.5); title('confounds','FontSize',10);
                                xlabel('Time in ms','FontSize',16)
                            end
                        end
                        
                        if strcmp(extra,'Adjusted')
                            if max(regressor) <= sum(LIMO.design.nb_conditions)
                                subplot(5,1,1:4);
                                plot(timevect,[ERP'-ERP2'],'LineWidth',1.5);
                                mytitle = sprintf('Adjusted ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                                assignin('base','Plotted_data',[ERP'-ERP2'])
                            else
                                subplot(5,1,1:4);
                                surf((squeeze(ERP(i,:,:)) - squeeze(tmp2(:,index)))');shading interp
                                mytitle = sprintf('Adjusted single trials \n sorted by regressor %g electrode %s (%g)', regressor, LIMO.data.chanlocs(electrode).labels, electrode);
                                assignin('base','Plotted_data',(squeeze(ERP(i,:,:)) - squeeze(tmp2(:,index)))')
                            end
                        elseif strcmp(extra,'Original')
                            if max(regressor) <= sum(LIMO.design.nb_conditions)
                                plot(timevect,ERP','LineWidth',1.5);
                                mytitle = sprintf('ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                                assignin('base','Plotted_data', ERP')
                            else
                                surf(squeeze(ERP(i,:,:))');shading interp
                                mytitle = sprintf('Original Single trials \n sorted by regressor %g electrode %s (%g)', regressor, LIMO.data.chanlocs(electrode).labels, electrode);
                                assignin('base','Plotted_data',squeeze(ERP(i,:,:))' )
                            end
                        else
                            if max(regressor) <= sum(LIMO.design.nb_conditions)
                                plot(timevect,ERP','LineWidth',1.5);
                                mytitle = sprintf('Modeled ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                                assignin('base','Plotted_data',ERP' )
                            else
                                subplot(1,1,1)
                                surf(squeeze(ERP(i,:,:))');shading interp
                                view(-19,38)
                                mytitle = sprintf('Modeled single trials \n sorted by regressor %g electrode %s (%g)', regressor, LIMO.data.chanlocs(electrode).labels, electrode);
                                assignin('base','Plotted_data', squeeze(ERP(i,:,:))')
                            end
                        end
                        
                        if max(regressor) <= sum(LIMO.design.nb_conditions)
                            grid on; axis tight; title(mytitle,'FontSize',19); drawnow;
                            v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                            set(gca,'FontSize',14);
                            ylabel('Amplitude in {\mu}V','FontSize',16)
                            if strcmp(extra,'Original')
                                xlabel('Time in ms','FontSize',16)
                            end
                        else
                            axis tight; title([mytitle,'   '],'FontSize',14); drawnow;
                            ylabel('Sorted trials','FontSize',16)
                            xlabel('Time in ms','FontSize',16)
                            zlabel('Amplitude in {\mu}V','FontSize',16)
                            labels_time_(LIMO, gca);
                            rotate3d
                        end
                    end
                end
            end
            
            
    end % closes switch
    
    
    
    
    
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % -------------------                               LEVEL 2
    % -------------------                            GROUP EFFECTs
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
elseif LIMO.Level == 2
    
    
    % ------------------------------
    %      Image and topoplot
    % ----------------------------
    if Type == 1 || Type == 2
        
        % special cases for ANCOVA
        % ---------------------------
        if strcmp(LIMO.design.name,'One-way ANCOVA one electrode') || strcmp(LIMO.design.name,'One-way ANCOVA all electrodes')
            % -----------------------------------------------------------------------------------------------------------------
            option = questdlg('which effect to plot','Plotting option','gp effect','covariates','gp effect');
            if strcmp(option,'gp effect')
                effect_nb = 1;
            else
                if size(N_way_covariance_analysis,3) == 2
                    effect_nb = 2;
                else
                    limo_review(LIMO);
                    effect_nb = eval(cell2mat(inputdlg('which covariate to plot','Plotting option')));
                    if length(effect_nb)>1
                        go =0;
                        while go==0
                            errordlg('please select only one regressor');
                            effect_nb = eval(cell2mat(inputdlg('which covariate to plot','Plotting option')))-(LIMO.design.nb_conditions-1);
                            if length(effect_nb) == 1 ; go =1; end
                        end
                    end
                    close('Review Design')
                end
            end
            [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,[],effect_nb);
            
            
        else  % all other designs
            % ----------------------
            [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,[],[]);
        end
        
        if isempty(M)
            errordlg('File not supported');
            return
        end
        
        if Type == 1
            %--------------------------
            % imagesc of the results
            %--------------------------
            if numel(size(M)) == 2
                if sum(mask(:)) == 0
                    warndlg('no values under threshold parameter','no significant effect');
                else
                    scale = M.*mask;
                    v = max(scale(:)); [e,f]=find(scale==v);
                    if min(scale(:))<0
                        scale(scale==0)=min(scale(:))+(min(scale(:))/10);
                    else
                        scale(scale==0)=NaN;
                    end
                    
                    figure; set(gcf,'Color','w');
                    ax(1) = subplot(2,3,[1 2 4 5]);
                    timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(M,2));
                    ratio = (LIMO.data.end*1000 - LIMO.data.start*1000) / size(M,2);
                    if LIMO.data.start < 0
                        frame_zeros = round(abs(LIMO.data.start*1000) / ratio);
                    end
                    imagesc(timevect,1:size(M,1),scale);
                    title([mytitle],'FontSize',18);
                    color_images_(scale,LIMO);
                    assignin('base','Plotted_data',scale)
                    
                    ax(2) = subplot(2,3,3);
                    topoplot(M(:,f),LIMO.data.chanlocs);
                    title(['topoplot @' num2str(timevect(f)) 'ms'],'FontSize',14)
                    
                    ax(3) = subplot(2,3,6);
                    plot(M(e,:),'LineWidth',3); grid on
                    mytitle = sprintf('time course @ \n electrode %s', LIMO.data.chanlocs(e).labels);
                    title(mytitle,'FontSize',14)
                    
                    update = 0;
                    while update ==0
                        try % use try so that if figure deleted no complain
                            [x,y,button]=ginput(1);
                        catch
                            update = 1; break
                        end
                        if button > 1
                            update = 1;
                        end
                        clickedAx = gca;
                        if clickedAx ~=ax(1)
                            disp('right click to exit')
                        else
                            subplot(2,3,3,'replace');
                            frame = round(x / ratio);
                            if frame < 0
                                frame = frame_zeros + frame;
                            end
                            topoplot(M(:,frame),LIMO.data.chanlocs);
                            title(['topoplot @ ' num2str(round(x)) 'ms'],'FontSize',14)
                            subplot(2,3,6,'replace'); y = round(y);
                            plot(M(y,:),'LineWidth',3); grid on
                            mytitle = sprintf('time course @ \n electrode %s', LIMO.data.chanlocs(y).labels);
                            title(mytitle,'FontSize',14)
                            clear x y button
                        end
                    end
                end
                
            else
                for i=1:size(M,3)-1
                    if sum(mask(:,:,i)) == 0
                        warndlg(['no values under threshold parameter',num2str(i)],'no significant effect');
                    else
                        figure('Name',['Parameter',num2str(i)]);set(gcf,'Color','w');
                        toplot = M(:,:,i);
                        scale = toplot.*mask(:,:,i);
                        if min(scale(:))<0
                            scale(scale==0)=min(scale(:))+(min(scale(:))/10);
                        else
                            scale(scale==0)=NaN;
                        end
                        imagesc(linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(toplot,2)),1:size(toplot,1),scale);
                        title([mytitle],'FontSize',20);
                        color_images_(scale,LIMO);
                        assignin('base','Plotted_data',scale)
                    end
                end
            end
            
            
        elseif Type == 2
            %--------------------------
            % topoplot
            %--------------------------
            if ~isempty(LIMO.design.electrode)  % not full scalp
                msgbox('Only one electrode found','No topoplot')
            elseif sum(mask(:)) == 0
                warndlg('no values under threshold','no significant effect');
            else
                EEG.data = M.*mask;
                EEG.setname = mytitle;
                EEG.pnts = size(EEG.data,2);
                EEG.xmin = LIMO.data.start;
                EEG.xmax = LIMO.data.end;
                EEG.times =  (LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000); % in sec;
                EEG.trials = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                EEG.nbchan = size(EEG.data,1);
                pop_topoplot(EEG);
                assignin('base','Plotted_data',EEG.data)
            end
        end
        
        
    elseif Type == 3
        
        %--------------------------
        % ERP
        %--------------------------
        
        if strcmp(LIMO.design.name,'one sample t-test one electrode') || strcmp(LIMO.design.name,'one sample t-test all electrodes')
            % --------------------------------------------------------------------------------------------------------------------------
            % one sample(electrodes,frames,[mean,std,n,t,p])
            % boot_one_sample(electrodes,frames,[sorted t, p],nboot)
            
            full_file = dir(FileName);
            load(full_file.name);
            [e,f,d]=size(one_sample);
            if e > 1
                electrode = inputdlg('which electrode to plot','Plotting option');
                if isempty(electrode) || strcmp(cell2mat(electrode),'')
                    [v,e] = max(one_sample(:,:,4)); [v,c]=max(v); electrode = e(c);
                else
                    electrode = eval(cell2mat(electrode));
                    if length(electrode) > 1
                        error('1 electrode only can be plotted')
                    elseif electrode > e
                        error('electrode number invalid')
                    end
                end
            else
                electrode = 1;
            end
            
            figure;set(gcf,'Color','w')
            timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
            plot(timevect,squeeze(one_sample(electrode,:,1)),'LineWidth',3);
            n = one_sample(electrode,:,3); s = one_sample(electrode,:,2);
            
            % check for boot parameters
            try
                MCC_data = sprintf('boot_%s',full_file.name); load(MCC_data);
                if MCC ==1  %% if not corrected but boot data are there ask
                    choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
                    if strcmp(choice,'use empirical p values')
                        robust =1;
                    else
                        robust = 0;
                    end
                else % MCC choosen and boot data are there
                    robust =1;
                end
            catch ME
                robust = 0;
            end
            
            if robust == 1
                lo = round(LIMO.design.nboot.*p./2);
                hi = LIMO.design.nboot - lo;
                sort_boot_one_sample = sort(squeeze(boot_one_sample(electrode,:,2,:)),2); % get T value under H0
                b              = squeeze(one_sample(electrode,:,1)) - (sort_boot_one_sample(:,lo+1)' .* (s ./ sqrt(n)));
                c              = squeeze(one_sample(electrode,:,1)) - (sort_boot_one_sample(:,hi)' .* (s ./ sqrt(n)));
                [M, mask, mytitle] = limo_stat_values(Type,full_file.name,p,MCC,LIMO,'use empirical p values',[]);
                if electrode == 1
                    if length(LIMO.design.electrode) == 1
                        mytitle = sprintf('Mean parameter value and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                    else
                        mytitle = sprintf('Mean parameter value and %g%% bootstrap CI \n optimized electrode', (1-p).*100);
                    end
                else
                    mytitle = sprintf('Mean parameter value and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                end
            else
                critical_value = tinv((1 - p / 2), n - 1) .* (one_sample(electrode,:,2) ./ sqrt(n)); % one_sample(:,:,2)=std
                b              = squeeze(one_sample(electrode,:,1)) - critical_value;
                c              = squeeze(one_sample(electrode,:,1)) + critical_value;
                [M, mask, mytitle] = limo_stat_values(Type,full_file.name,p,MCC,LIMO,'use theoretical p values',[]);
                if electrode == 1
                    if length(LIMO.design.electrode) == 1
                        mytitle = sprintf('Mean parameter value and %g%% CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                    else
                        mytitle = sprintf('Mean parameter value and %g%% CI \n optimized electrode', (1-p).*100);
                    end
                else
                    mytitle = sprintf('Mean parameter value and %g%% CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                end
            end
            
            hold on; plot(timevect,b,'--r','LineWidth',2);plot(timevect,c,'--r','LineWidth',2)
            xlabel('Time in ms','FontSize',16)
            ylabel('Amplitude in {\mu}V per std','FontSize',16)
            grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
            v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
            set(gca,'FontSize',14);
            sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
            plot(timevect,sig.*v(3),'r.','LineWidth',3);
            
            plotted_data = NaN(size(one_sample,2),4);
            plotted_data(:,1) = c; % upper CI
            plotted_data(:,2) = squeeze(one_sample(electrode,:,1)); % data
            plotted_data(:,3) = b; % lower CI
            plotted_data(:,4) = sig; % Sig
            assignin('base','Plotted_data',plotted_data)
            
            
        elseif strcmp(LIMO.design.name,'two samples t-test one electrode') || strcmp(LIMO.design.name,'two samples t-test all electrodes') ...
                % --------------------------------------------------------------------------------------------------------------------------
            
            % two samples(electrodes, frames, [mean diff, dfe, std gp1, std gp2, t, p])
            % boot_two_samples(electrodes,frames,[sorted differences, sorted t, p],nboot)
            
            full_file = dir('two_samples*');
            load(full_file.name);
            [e,f,d]=size(two_samples_ttest);
            if e > 1
                electrode = inputdlg('which electrode to plot','Plotting option');
                if isempty(electrode) || strcmp(cell2mat(electrode),'')
                    [v,e] = max(two_samples_ttest(:,:,4)); [v,c]=max(v); electrode = e(c);
                else
                    electrode = eval(cell2mat(electrode));
                    if length(electrode) > 1
                        error('1 electrode only can be plotted')
                    elseif electrode > e
                        error('electrode number invalid')
                    end
                end
            else
                electrode = 1;
            end
            
            figure;set(gcf,'Color','w')
            timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
            plot(timevect,squeeze(two_samples_ttest(electrode,:,1)),'LineWidth',3);
            
            % check for boot parameters
            try
                MCC_data = sprintf('boot_%s',full_file.name); load(MCC_data);
                if MCC ==1  %% if not corrected but boot data are there ask
                    choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
                    if strcmp(choice,'use empirical p values')
                        robust =1;
                    else
                        robust = 0;
                    end
                else % MCC choosen and boot data are there
                    robust =1;
                end
            catch ME
                robust = 0;
            end
            
            if robust == 1
                
                % percentile bootstrap 95% CI of the mean difference
                lo = round(LIMO.design.nboot.*p./2);
                hi = LIMO.design.nboot - lo;
                sort_boot_two_samples_ttest = sort(boot_two_samples_ttest,4); % sort boot
                b = sort_boot_two_samples_ttest(electrode,:,1,lo+1); % take T
                c = sort_boot_two_samples_ttest(electrode,:,1,hi);
                [M, mask, mytitle] = limo_stat_values(Type,full_file.name,p,MCC,LIMO,'use empirical p values',[]); % 'use empirical p values'
                if electrode == 1
                    if length(LIMO.design.electrode) == 1
                        mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                    else
                        mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI \n optimized electrode', (1-p).*100);
                    end
                else
                    mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                end
            else
                se = two_samples_ttest(electrode,:,5) ./ two_samples_ttest(electrode,:,1); % t = m ./ se
                dfe = two_samples_ttest(electrode,:,2);
                critical_value = tinv(1 - p ./ 2, dfe) .* se;
                b = two_samples_ttest(electrode,:,1) - critical_value;
                c = two_samples_ttest(electrode,:,1) + critical_value;
                [M, mask, mytitle] = find_values(Type,full_file.name,p,MCC,LIMO,'use theoretical p values',[]);
                if electrode == 1
                    if length(LIMO.design.electrode) == 1
                        mytitle = sprintf('Mean parameter difference between groups and %g%% CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                    else
                        mytitle = sprintf('Mean parameter difference between groups and %g%% CI \n optimized electrode', (1-p).*100);
                    end
                else
                    mytitle = sprintf('Mean parameter difference between groups and %g%% CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                end
            end
            
            hold on; plot(timevect,b,'--r','LineWidth',2);plot(timevect,c,'--r','LineWidth',2)
            xlabel('Time in ms','FontSize',16)
            ylabel('Amplitude in {\mu}V per std','FontSize',16)
            grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
            v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
            set(gca,'FontSize',14);
            sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
            plot(timevect,sig.*v(3),'r.','LineWidth',3);
            
            plotted_data = NaN(size(two_samples_ttest,2),4);
            plotted_data(:,1) = c; % upper CI
            plotted_data(:,2) = squeeze(two_samples_ttest(electrode,:,1)); % data
            plotted_data(:,3) = b; % lower CI
            plotted_data(:,4) = sig; % Sig
            assignin('base','Plotted_data',plotted_data)
            
            
        elseif strcmp(LIMO.design.name,'paired t-test one electrode') || strcmp(LIMO.design.name,'paired t-test all electrodes')
            % --------------------------------------------------------------------------------------------------------------------------
            
            % works exactly like the one sample t-test because the test is a one
            % sample on the difference between pairs
            
            load('paired_ttest');
            [e,f,d]=size(paired_ttest);
            if e > 1
                electrode = inputdlg('which electrode to plot','Plotting option');
                if isempty(electrode) || strcmp(cell2mat(electrode),'')
                    [v,e] = max(paired_ttest(:,:,4)); [v,c]=max(v); electrode = e(c);
                else
                    electrode = eval(cell2mat(electrode));
                    if length(electrode) > 1
                        error('1 electrode only can be plotted')
                    elseif electrode > e
                        error('electrode number invalid')
                    end
                end
            else
                electrode = 1;
            end
            
            figure;set(gcf,'Color','w')
            timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in ms
            plot(timevect,squeeze(paired_ttest(electrode,:,1)),'LineWidth',3);
            n = paired_ttest(electrode,:,3);
            s = paired_ttest(electrode,:,2); % std
            
            
            % check for boot parameters
            try
                load('boot_paired_ttest')
                if MCC ==1  %% if not corrected but boot data are there ask
                    choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
                    if strcmp(choice,'use empirical p values')
                        robust = 1;
                    else
                        robust = 0;
                    end
                else % MCC choosen and boot data are there
                    robust =1;
                end
            catch ME
                robust = 0;
            end
            
            if robust == 1
                lo = round(LIMO.design.nboot.*(p/2));
                hi = LIMO.design.nboot - lo;
                Tsorted  = sort(squeeze(boot_paired_ttest(electrode,:,1,:)),2); % take the mean diff. under H1 along bootstrap dimension
                b = Tsorted(:,lo+1);
                c = Tsorted(:,hi);
                [M, mask, mytitle] = limo_stat_values(Type,'paired_ttest',p,MCC,LIMO,'use empirical p values',[]);
                if electrode == 1
                    if length(LIMO.design.electrode) == 1
                        mytitle = sprintf('Mean parameter difference \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                    else
                        mytitle = sprintf('Mean parameter difference \n and %g%% bootstrap CI \n optimized electrode', (1-p).*100);
                    end
                else
                    mytitle = sprintf('Mean parameter difference \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                end
            else
                critical_value = tinv((1 - p / 2), n - 1) .* (paired_ttest(electrode,:,2) ./ sqrt(n)); % paired_ttest(electrode,:,2) = std
                b              = squeeze(paired_ttest(electrode,:,1)) - critical_value;
                c              = squeeze(paired_ttest(electrode,:,1)) + critical_value;
                [M, mask, mytitle] = limo_stat_values(Type,'paired_ttest',p,MCC,LIMO,'use theoretical p values',[]);
                if electrode == 1
                    if length(LIMO.design.electrode) == 1
                        mytitle = sprintf('Mean parameter difference and %g%% CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(1).labels,LIMO.design.electrode);
                    else
                        mytitle = sprintf('Mean parameter difference value %g%% CI \n optimized electrode', (1-p).*100);
                    end
                else
                    mytitle = sprintf('Mean parameter difference and %g%% CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                end
            end
            
            hold on; plot(timevect,b,'--r','LineWidth',2);plot(timevect,c,'--r','LineWidth',2)
            xlabel('Time in ms','FontSize',16)
            ylabel('Amplitude in {\mu}V per std','FontSize',16)
            grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
            v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
            set(gca,'FontSize',14);
            sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
            plot(timevect,sig.*v(3),'r.','LineWidth',3);
            
            plotted_data = NaN(size(paired_ttest,2),4);
            plotted_data(:,1) = c; % upper CI
            plotted_data(:,2) = squeeze(paired_ttest(electrode,:,1)); % data
            plotted_data(:,3) = b; % lower CI
            plotted_data(:,4) = sig; % Sig
            assignin('base','Plotted_data',plotted_data)
            
            
        elseif strcmp(LIMO.design.name,'regression analysis one electrode') || strcmp(LIMO.design.name,'regression analysis all electrodes')
            % --------------------------------------------------------------------------------------------------------------------------
            
            load('Betas.mat');
            FileName = dir('regression*');
            FileName = FileName.name;
            load(FileName);
            [e,f,d]=size(Betas);
            
            % check with regressor to plot
            if d > 2
                limo_review(LIMO);
                effect_nb = eval(cell2mat(inputdlg('which regressor to plot','Plotting option')));
                if length(effect_nb)>1
                    go =0;
                    while go==0
                        errordlg('please select only one regressor');
                        effect_nb = eval(cell2mat(inputdlg('which covariate to plot','Plotting option')));
                        if length(effect_nb) == 1 ; go =1; end
                    end
                end
                close('Review Design')
            else
                effect_nb = 1;
            end
            
            % check which electrode to plot
            % regression = (electrodes, frames, nb_regressors, [F / p])
            if e > 1
                electrode = inputdlg('which electrode to plot','Plotting option');
                if isempty(electrode) || strcmp(cell2mat(electrode),'')
                    clear electrode; [v,e] = max(regression(:,:,effect_nb,1)); [v,c]=max(v); electrode = e(c);
                else
                    electrode = eval(cell2mat(electrode));
                    if length(electrode) > 1
                        error('1 electrode only can be plotted')
                    elseif electrode > e
                        error('electrode number invalid')
                    end
                end
            else
                if length(LIMO.design.electrode) == 1
                    electrode = LIMO.design.electrode;
                else
                    electrode = 1;  % accomodates the fact that all matrices have the electrode dim (even = 1)
                end
            end
            
            
            % check for boot parameters
            try
                FileName2 = dir('boot_H1*');
                FileName2 = FileName2.name;
                load(FileName2);
                if MCC ==1  %% if not corrected but boot data are there ask
                    choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
                    if strcmp(choice,'use empirical p values')
                        robust =1;
                    else
                        robust = 0;
                    end
                else % MCC choosen and boot data arte there
                    robust =1;
                end
            catch ME
                robust = 0;
            end
            
            
            % get the CI
            data = NaN(size(Betas,1),size(Betas,2),3);
            if robust == 1
                data(:,:,2) = Betas(:,:,effect_nb);
                boot_data = squeeze(boot_H1_Betas(:,:,effect_nb,:));
                if numel(size(boot_data)) == 2  % if one electrode only
                    boot_data = ones(1,size(boot_data,1),size(boot_data,2));
                    boot_data(1,:,:) = squeeze(boot_H1_Betas(:,:,effect_nb,:));
                end
                boot_data = sort(boot_data,3);
                if size(LIMO.design.X,2) == 2 && LIMO.design.nboot == 599  % when there is only 1 regressor (dim = 2 because of the intercept)
                    % special adjustments for the slope estimate based on 600 bootstraps, based on Wilcox, 2005, p. 418-419.
                    if size(LIMO.design.X,1) < 40
                        a=6; c=593;
                    elseif size(LIMO.design.X,1) >=40 && size(LIMO.design.X,1) < 80
                        a=7; c=592;
                    elseif size(LIMO.design.X,1) >=80 && size(LIMO.design.X,1) < 180
                        a=10; c=589;
                    elseif size(LIMO.design.X,1) >=180 && size(LIMO.design.X,1) < 250
                        a=13; c=586;
                    else
                        a=15; c=584;
                    end
                    data(:,:,1) = boot_data(:,:,a+1);
                    data(:,:,3) = boot_data(:,:,c);
                else % we have several regressors - use Bonferroni approach for simultaneous prob coverage
                    low = round((p*LIMO.design.nboot) / (2*size(LIMO.design.X,2)));
                    high = LIMO.design.nboot - low;
                    data(:,:,1) = boot_data(:,:,low);
                    data(:,:,3) = boot_data(:,:,high);
                end
                
                % figure
                figure;set(gcf,'Color','w')
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                plot(timevect,squeeze(data(electrode,:,2)),'LineWidth',3);
                hold on; plot(timevect,squeeze(data(electrode,:,1)),'--r','LineWidth',2);
                plot(timevect,squeeze(data(electrode,:,3)),'--r','LineWidth',2)
                xlabel('Time in ms','FontSize',16)
                ylabel('Amplitude (arbitrary unit)','FontSize',16)
                [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,'use empirical p values',effect_nb); % note that if MCC ~=1 we don't care about the choice option
                if electrode ==1 && length(LIMO.design.electrode) ~=1 % ie this is an optimized electrode
                    mytitle = {['Regression coefficient ' num2str(effect_nb) ' and ' num2str((1-p).*100) '% bootstrap CI'];['optimized electrode']};
                else
                    mytitle = {['Regression coefficient ' num2str(effect_nb) ' and ' num2str((1-p).*100) '% bootstrap CI'];['at electrode ' LIMO.data.chanlocs(electrode).labels '(' num2str(electrode) ')']};
                end
                grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                plot(timevect,sig.*v(3),'r.','LineWidth',3);
                
                plotted_data = NaN(size(data,2),4);
                plotted_data(:,1) = squeeze(data(electrode,:,1)); % upper CI
                plotted_data(:,2) = squeeze(data(electrode,:,2)); % data
                plotted_data(:,3) = squeeze(data(electrode,:,3)); % lower CI
                plotted_data(:,4) = sig; % Sig
                assignin('base','Plotted_data',squeeze(plotted_data));  % with the squeeze if N=1, reduces the matrix
                
            else
                load Yr; load Res;
                df = sum((Yr ~= NaN),3)-2; t = tcdf(1-p,df);
                sigma2 = nansum((Res.^2 ./ repmat(df,[1 1 size(Res,3)])),3); clear Yr Res
                X = zscore(LIMO.design.X(:,effect_nb));
                data(:,:,2) = Betas(:,:,effect_nb);
                data(:,:,1) = data(:,:,2) - t.*sqrt(sigma2 ./ norm(X).^2);
                data(:,:,3) = data(:,:,2) + t.*sqrt(sigma2 ./ norm(X).^2);
                
                % figure
                figure;set(gcf,'Color','w')
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                plot(timevect,data(electrode,:,2),'LineWidth',3);
                hold on; plot(timevect,data(electrode,:,1),'--r','LineWidth',2);
                plot(timevect,data(electrode,:,3),'--r','LineWidth',2)
                xlabel('Time in ms','FontSize',16)
                ylabel('Amplitude (arbitrary unit)','FontSize',16)
                [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,'use theoretical p values',effect_nb);
                if electrode ==1 && length(LIMO.design.electrode) ~=1 % ie this is an optimized electrode
                    mytitle = {['Regression coefficient ' num2str(effect_nb) ' and ' num2str((1-p).*100) '% bootstrap CI'];['optimized electrode']};
                else
                    mytitle = {['Regression coefficient ' num2str(effect_nb) ' and ' num2str((1-p).*100) '% bootstrap CI'];['at electrode ' LIMO.data.chanlocs(electrode).labels '(' num2str(electrode) ')']};
                end
                grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                plot(timevect,sig.*v(3),'r.','LineWidth',3);
                
                plotted_data = NaN(size(data,2),4);
                plotted_data(:,1) = squeeze(data(electrode,:,1)); % upper CI
                plotted_data(:,2) = squeeze(data(electrode,:,2)); % data
                plotted_data(:,3) = squeeze(data(electrode,:,3)); % lower CI
                plotted_data(:,4) = sig; % Sig
                assignin('base','Plotted_data',squeeze(plotted_data));
            end
            
            
        elseif strcmp(LIMO.design.name,'Repeated measure ANOVA all electrodes') || strcmp(LIMO.design.name,'Repeated measure ANOVA one electrode') || ...
                strcmp(LIMO.design.name,'Mixed ANOVA all electrodes') || strcmp(LIMO.design.name,'Mixed ANOVA ANOVA one electrode')
            % --------------------------------------------------------------------------------------------------------------------------------------
            
            
            if strcmp(LIMO.design.name,'Mixed ANOVA all electrodes') || strcmp(LIMO.design.name,'Mixed ANOVA ANOVA one electrode')
                type = questdlg('Several effects to be visualized','plotting option','Repeated measures','Group','Group*Repeated measures','Gp*Repeated measures');
            else
                type ='Repeated measures';
            end
            
            % --------------------------------
            if strcmp(type,'Repeated measures')
                % ------------------------------
                load('Rep_ANOVA');
                [e,f,n,d]=size(Rep_ANOVA);
                
                if numel(size(Rep_ANOVA)) > 3
                    effect_nb = eval(cell2mat(inputdlg({['which factor to plot max ' num2str(n)]},'Plotting option')));
                    if effect_nb > n % nb of n parameters
                        go =0;
                        while go==0
                            errordlg({['there is only ' num2str(c) ' factors']});
                            effect_nb = eval(cell2mat(inputdlg('which factor to plot','Plotting option')));
                            if effect_nb <= n  ; go =1; end
                        end
                    end
                    Rep_ANOVA = squeeze(Rep_ANOVA(:,:,effect_nb,:));
                else
                    effect_nb = 1;
                end
                
                if e > 1
                    electrode = inputdlg('which electrode to plot','Plotting option');
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        [v,e] = max(Rep_ANOVA(:,:,1)); [v,c]=max(v); electrode = e(c);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > e
                            error('electrode number invalid')
                        end
                    end
                else
                    electrode = 1;
                end
                
                figure;set(gcf,'Color','w')
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in ms
                % the data to plot are the difference in Yr given
                % LIMO.design.C (see limo_rep_anova)
                C = LIMO.design.C{effect_nb};
                load Yr
                Data = squeeze(Yr(electrode,:,:,:));
                % compute differences between pairs using C
                for time=1:size(Data,1)
                    avg(:,time) = mean(C*squeeze(Data(time,:,:))',2);
                end
                plot(timevect,avg,'LineWidth',3);
                
                % check for boot parameters
                try
                    load('boot_index')
                    if MCC ==1  %% if not corrected but boot data are there ask
                        choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
                        if strcmp(choice,'use empirical p values')
                            robust = 1;
                        else
                            robust = 0;
                        end
                    else % MCC choosen and boot data are there
                        robust =1;
                    end
                catch ME
                    robust = 0;
                end
                
                if robust == 1
                    lo = round(LIMO.design.nboot.*(p/2));
                    hi = LIMO.design.nboot - lo;
                    % compute the mean difference under H1
                    for boot = 1:size(boot_index,2)
                        Data = squeeze(Yr(electrode,:,boot_index(:,boot),:));
                        for time=1:size(Data,1)
                            avgb(:,time,boot) = mean(C*squeeze(Data(time,:,:,:))',2);
                        end
                    end
                    Meansorted  = sort(avgb,numel(size(avgb))); % take the mean diff. along bootstrap dimension
                    [M, mask, mytitle] = limo_stat_values(Type,'Rep_ANOVA',p,MCC,LIMO,'use empirical p values',effect_nb); % no need to use effect_nb since we squeeze(Rep_ANOVA) above
                    sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                    
                    if numel(size(avgb)) == 2
                        b = Meansorted(:,lo+1);
                        c = Meansorted(:,hi);
                        if electrode == 1
                            if length(LIMO.design.electrode) == 1
                                mytitle = sprintf('Mean parameters difference \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                            else
                                mytitle = sprintf('Mean parameters difference \n and %g%% bootstrap CI \n optimized electrode', (1-p).*100);
                            end
                        else
                            mytitle = sprintf('Mean parameters difference \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                        end
                        hold on; plot(timevect,b,'--r','LineWidth',2);plot(timevect,c,'--r','LineWidth',2)
                        ylabel('Amplitude of effect (arbitrary unit)','FontSize',16)
                        plotted_data = NaN(size(avg,2),4);
                        plotted_data(:,1) = c; % upper CI
                        plotted_data(:,2) = avg; % data
                        plotted_data(:,3) = b; % lower CI
                        plotted_data(:,4) = sig; % Sig
                        
                    else
                        b = Meansorted(:,:,lo+1);
                        c = Meansorted(:,:,hi);
                        if electrode == 1
                            if length(LIMO.design.electrode) == 1
                                mytitle = sprintf('Mean parameter differences \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                            else
                                mytitle = sprintf('Mean parameter differences \n and %g%% bootstrap CI \n optimized electrode', (1-p).*100);
                            end
                        else
                            mytitle = sprintf('Mean parameter differences \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                        end
                        hold on; plot(timevect,b,'--','LineWidth',2);plot(timevect,c,'--','LineWidth',2)
                        ylabel('Amplitude of effects (arbitrary unit)','FontSize',16)
                        plotted_data = NaN(size(avg,1),size(avg,2),4);
                        plotted_data(:,:,1) = c; % upper CI
                        plotted_data(:,:,2) = avg; % data
                        plotted_data(:,:,3) = b; % lower CI
                        plotted_data(:,:,4) = repmat(sig,size(avg,1),1); % Sig
                    end
                    
                    xlabel('Time in ms','FontSize',16)
                    grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                    v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                    set(gca,'FontSize',14);
                    plot(timevect,sig.*v(3),'r.','LineWidth',3);
                    assignin('base','Plotted_data',plotted_data)
                    
                else
                    % CI
                    S = NaN(size(Data,1),size(Data,3),size(Data,3));
                    for time = 1:size(Data,1)
                        S(time,:,:) = cov(squeeze(Data(time,:,:))); % covariance to account for spericity
                    end
                    k = size(C,1) + 1;
                    alpha = (2*p) / (k*(k-1)); % = p if k = 2
                    df = size(Data,3)-1;
                    dfe = size(Data,2)-size(Data,3)+1;
                    c = avg + 2*finv(alpha,df,dfe).*(sqrt(C*squeeze(S(time,:,:))*C')); % uses Bonferoni inequality
                    b = avg - 2*finv(alpha,df,dfe).*(sqrt(C*squeeze(S(time,:,:))*C'));
                    % Sig
                    [M, mask, mytitle] = limo_stat_values(Type,'Rep_ANOVA',p,MCC,LIMO,'use theoretical p values',effect_nb); % no need to use effect_nb since we squeeze(Rep_ANOVA) above
                    sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                    v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                    hold on; plot(timevect,sig.*v(3),'r.','LineWidth',3); drawnow;
                    plot(timevect,b','--r','LineWidth',2);plot(timevect,c','--r','LineWidth',2)
                    mytitle = sprintf('Mean parameters difference\n at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels,electrode);
                    xlabel('Time in ms','FontSize',16)
                    ylabel('Amplitude in {\mu}V per std','FontSize',16)
                    grid on; axis tight; title(mytitle,'FontSize',20);
                    set(gca,'FontSize',14);
                    
                    plotted_data = NaN(size(avg,2),4);
                    plotted_data(:,1) = c; % upper CI
                    plotted_data(:,2) = avg; % data
                    plotted_data(:,3) = b; % lower CI
                    plotted_data(:,4) = sig; % Sig
                    assignin('base','Plotted_data',plotted_data)
                    
                end
                
                % ----------------------
            elseif strcmp(type,'Group')
                % ----------------------
                load('Rep_ANOVA_Gp_effect')
                [e,f,d]=size(Rep_ANOVA_Gp_effect);
                
                % check electrode to plot
                if e > 1
                    electrode = inputdlg('which electrode to plot','Plotting option');
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        clear electrode; [v,e] = max(Rep_ANOVA_Gp_effect(:,:,1)); [v,c]=max(v); electrode = e(c);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > e
                            error('electrode number invalid')
                        end
                    end
                else
                    if length(LIMO.design.electrode) == 1
                        electrode = LIMO.design.electrode;
                    else
                        electrode = 1;  % accomodates the fact that all matrices have the electrode dim (even = 1)
                    end
                end
                
                % compute the pair-wise differences and plot
                figure;set(gcf,'Color','w')
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                load Yr; Data = mean(squeeze(Yr(electrode,:,:,:)),3);
                combinations = nchoosek(1:max(LIMO.design.gp),2);
                for d = 1:size(combinations,1)
                    Effect(:,d) = mean(Data(:,find(LIMO.design.gp == combinations(d,1))),2) - mean(Data(:,find(LIMO.design.gp == combinations(d,2))),2);
                end
                plot(timevect,Effect,'LineWidth',3);
                
                
                % check for boot parameters
                try
                    load('boot_index')
                    if MCC ==1  %% if not corrected but boot data are there ask
                        choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
                        if strcmp(choice,'use empirical p values')
                            robust = 1;
                        else
                            robust = 0;
                        end
                    else % MCC choosen and boot data are there
                        robust =1; choice = 'use empirical p values';
                    end
                catch ME
                    robust = 0; choice = 'use theoretical p values';
                end
                [M, mask, mytitle] = limo_stat_values(Type,'Rep_ANOVA_Gp_effect',p,MCC,LIMO,choice,[]);
                
                % CI
                if robust == 1
                    lo = round(LIMO.design.nboot.*p./2);
                    hi = LIMO.design.nboot - lo;
                    for boot = 1:LIMO.design.nboot
                        for d = 1:size(combinations,1)
                            Data = mean(squeeze(Yr(electrode,:,boot_index(:,boot),:)),3);
                            Effectb(:,d,boot) = mean(Data(:,find(LIMO.design.gp == combinations(d,1))),2) - mean(Data(:,find(LIMO.design.gp == combinations(d,2))),2);
                        end
                    end
                    sort_boot_H1 = sort(Effectb,3); % sort boot
                    b = sort_boot_H1(:,:,lo+1);
                    c = sort_boot_H1(:,:,hi);
                    if electrode == 1
                        if length(LIMO.design.electrode) == 1
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI optimized electrode', (1-p).*100);
                        end
                    else
                        mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                    end
                    
                else
                    
                    Y = mean(squeeze(Yr(electrode,:,:,:)),3);
                    N = length(LIMO.design.gp);
                    X = zeros(size(Y,2),LIMO.design.nb_conditions+1);
                    X(:,end) = 1; index = 1;
                    for i=1:LIMO.design.nb_conditions
                        X(find(LIMO.design.gp == i),i) = 1;
                    end
                    Res = (Y' - X*(pinv(X)*Y'))';
                    df = sum((Y ~= NaN),2)-rank(X); t = tcdf(1-p,df);
                    sigma2 = nansum((Res.^2 ./ repmat(df,[1 size(Res,2)])),2);
                    v = t.*sqrt(sigma2 ./ norm(X(:,1:end-1)).^2);
                    b = Effect - v(electrode,:);
                    c = Effect + v(electrode,:);
                    if electrode == 1
                        if length(LIMO.design.electrode) == 1
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% CI optimized electrode', (1-p).*100);
                        end
                    else
                        mytitle = sprintf('Mean parameter difference between groups \n and %g%% CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                    end
                end
                
                hold on; plot(timevect,b,'--r','LineWidth',3);
                plot(timevect,c,'--r','LineWidth',3);
                xlabel('Time in ms','FontSize',16)
                ylabel('Amplitude (arbitrary unit)','FontSize',16)
                grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                hold on; plot(timevect,sig.*v(3),'r.','LineWidth',3);
                
                plotted_data(:,:,1) = c;
                plotted_data(:,:,2) = Effect;
                plotted_data(:,:,3) = b;
                plotted_data(:,:,4) = repmat(sig,size(Effect,2),1);
                assignin('base','Plotted_data',squeeze(plotted_data));
                
                
                
                % -------------------------
            else % Gp * Repeated measures
                % ------------------------
                
                load('Rep_ANOVA_Interaction_with_gp');
                Rep_ANOVA = Rep_ANOVA_Interaction_with_gp; % just to use a shorter name
                clear Rep_ANOVA_Interaction_with_gp;
                [e,f,n,d]=size(Rep_ANOVA);
                
                if numel(size(Rep_ANOVA)) > 3
                    effect_nb = eval(cell2mat(inputdlg({['which factor to plot max ' num2str(n)]},'Plotting option')));
                    if effect_nb > n % nb of n parameters
                        go =0;
                        while go==0
                            errordlg({['there is only ' num2str(c) ' factors']});
                            effect_nb = eval(cell2mat(inputdlg('which factor to plot','Plotting option')));
                            if effect_nb <= n  ; go =1; end
                        end
                    end
                    Rep_ANOVA = squeeze(Rep_ANOVA(:,:,effect_nb,:));
                else
                    effect_nb = 1;
                end
                
                if e > 1
                    electrode = inputdlg('which electrode to plot','Plotting option');
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        [v,e] = max(Rep_ANOVA(:,:,1)); [v,c]=max(v); electrode = e(c);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > e
                            error('electrode number invalid')
                        end
                    end
                else
                    electrode = 1;
                end
                
                figure;set(gcf,'Color','w')
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in ms
                % the data to plot are the difference in Yr given
                % LIMO.design.C (see limo_rep_anova)
                C = LIMO.design.C{effect_nb};
                load Yr; Data = squeeze(Yr(electrode,:,:,:));
                % compute differences between pairs using C but for each gp
                for g=1:LIMO.design.nb_conditions
                    for time=1:size(Data,1)
                        avg(:,time,g) = mean(C*squeeze(Data(time,find(LIMO.design.gp == g),:))',2);
                    end
                end
                
                hold on;
                for w=1:size(avg,1)
                    plot(timevect,squeeze(avg(w,:,:)),'LineWidth',3);
                end
                
                % check for boot parameters
                try
                    load('boot_index')
                    if MCC ==1  %% if not corrected but boot data are there ask
                        choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
                        if strcmp(choice,'use empirical p values')
                            robust = 1;
                        else
                            robust = 0;
                        end
                    else % MCC choosen and boot data are there
                        robust =1;
                    end
                catch ME
                    robust = 0;
                end
                
                if robust == 1
                    lo = round(LIMO.design.nboot.*(p/2));
                    hi = LIMO.design.nboot - lo;
                    % compute the mean difference under H1
                    for boot = 1:size(boot_index,2)
                        Data = squeeze(Yr(electrode,:,boot_index(:,boot),:));
                        for g=1:LIMO.design.nb_conditions
                            for time=1:size(Data,1)
                                avgb(:,time,g,boot) = mean(C*squeeze(Data(time,find(LIMO.design.gp==g),:))',2);
                            end
                        end
                    end
                    Meansorted  = sort(avgb,numel(size(avgb))); % take the mean diff. along bootstrap dimension
                    b = squeeze(Meansorted(:,:,:,lo+1));
                    c = squeeze(Meansorted(:,:,:,hi));
                    if electrode == 1
                        if length(LIMO.design.electrode) == 1
                            mytitle = sprintf('Mean parameter differences \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Mean parameter differences \n and %g%% bootstrap CI \n optimized electrode', (1-p).*100);
                        end
                    else
                        mytitle = sprintf('Mean parameter differences \n and %g%% bootstrap CI \n at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                    end
                    [M, mask, mytitle] = limo_stat_values(Type,'Rep_ANOVA',p,MCC,LIMO,'use empirical p values',effect_nb); % no need to use effect_nb since we squeeze(Rep_ANOVA) above
                    hold on;
                    plot(timevect,b,'--','LineWidth',2);plot(timevect,c,'--','LineWidth',2)
                    ylabel('Amplitude of effects (arbitrary unit)','FontSize',16)
                    xlabel('Time in ms','FontSize',16)
                    grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                    v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                    set(gca,'FontSize',14);
                    sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                    plot(timevect,sig.*v(3),'r.','LineWidth',3);
                    
                    plotted_data = NaN(size(avg,1),size(avg,2),size(avg,3),4);
                    plotted_data(:,:,:,1) = c; % upper CI
                    plotted_data(:,:,:,2) = avg; % data
                    plotted_data(:,:,:,3) = b; % lower CI
                    plotted_data(:,:,:,4) = repmat(sig,size(avg,3),1)'; % Sig
                    assignin('base','Plotted_data',plotted_data)
                    
                else
                    % CI
                    k = size(C,1) + 1;
                    alpha = (2*p) / (k*(k-1)); % = p if k = 2
                    df = size(Data,3)-1;
                    dfe = size(Data,2)-size(Data,3)+1;
                    S = NaN(size(Data,1),size(Data,3),size(Data,3));
                    for g=1:LIMO.design.nb_conditions
                        for time = 1:size(Data,1)
                            S(time,:,:) = cov(squeeze(Data(time,find(LIMO.design.gp == g),:))); % covariance to account for spericity
                        end
                        c{g} = avg(:,:,g) + 2*finv(alpha,df,dfe).*(sqrt(C*squeeze(S(time,:,:))*C')); % uses Bonferoni inequality
                        b{g} = avg(:,:,g) - 2*finv(alpha,df,dfe).*(sqrt(C*squeeze(S(time,:,:))*C'));
                    end
                    % Sig
                    [M, mask, mytitle] = limo_stat_values(Type,'Rep_ANOVA',p,MCC,LIMO,'use theoretical p values',effect_nb); % no need to use effect_nb since we squeeze(Rep_ANOVA) above
                    sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                    v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                    hold on; plot(timevect,sig.*v(3),'r.','LineWidth',3); drawnow;
                    plot(timevect,cell2mat(b'),'--','LineWidth',2);
                    plot(timevect,cell2mat(c')','--','LineWidth',2)
                    mytitle = sprintf('Mean parameters difference\n at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels,electrode);
                    xlabel('Time in ms','FontSize',16)
                    ylabel('Amplitude in {\mu}V per std','FontSize',16)
                    grid on; axis tight; title(mytitle,'FontSize',20);
                    set(gca,'FontSize',14);
                    
                    plotted_data = NaN(size(avg,1),size(avg,2),size(avg,3),4);
                    plotted_data(:,:,:,1) = c; % upper CI
                    plotted_data(:,:,:,2) = avg; % data
                    plotted_data(:,:,:,3) = b; % lower CI
                    plotted_data(:,:,:,4) = repmat(sig,size(avg,3),1)'; % Sig
                    assignin('base','Plotted_data',plotted_data)
                    
                end
            end
            
            
        elseif strcmp(LIMO.design.name,'One-way ANOVA one electrode') || strcmp(LIMO.design.name,'One-way ANOVA all electrodes')
            % ------------------------------------------------------------------------------------------------------------------
            
            load('Betas.mat');
            load('One_way_variance_analysis')
            
            [e,f,d]=size(Betas);
            if e > 1
                electrode = inputdlg('which electrode to plot','Plotting option');
                if isempty(electrode) || strcmp(cell2mat(electrode),'')
                    clear electrode; [v,e] = max(N_way_variance_analysis(:,:,1)); [v,c]=max(v); electrode = e(c);
                else
                    electrode = eval(cell2mat(electrode));
                    if length(electrode) > 1
                        error('1 electrode only can be plotted')
                    elseif electrode > e
                        error('electrode number invalid')
                    end
                end
            else
                if length(LIMO.design.electrode) == 1
                    electrode = LIMO.design.electrode;
                else
                    electrode = 1;  % accomodates the fact that all matrices have the electrode dim (even = 1)
                end
            end
            
            
            % check for boot parameters
            try
                load('boot_H1_Betas')
                if MCC ==1  %% if not corrected but boot data are there ask
                    choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
                    if strcmp(choice,'use empirical p values')
                        robust = 1;
                    else
                        robust = 0;
                    end
                else % MCC choosen and boot data are there
                    robust =1; choice = 'use empirical p values';
                end
            catch ME
                robust = 0;choice = 'use theoretical p values';
            end
            
            type = questdlg('2 types of plots are available','plotting option','plot all parameters','plot the difference between all parameters','plot the difference between all parameters');
            
            
            if strcmp(type,'plot the difference between all parameters')
                
                figure;set(gcf,'Color','w')
                [M, mask, mytitle] = limo_stat_values(Type,'One_way_variance_analysis',p,MCC,LIMO,choice,[]);
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                c = [eye(LIMO.design.nb_conditions-1) ones(LIMO.design.nb_conditions-1,1).*-1];
                C = limo_contrast_checking(LIMO.dir, LIMO.design.X, c);
                Effect = C*squeeze(Betas(electrode,:,:))';
                
                if robust == 1
                    H1 = squeeze(boot_H1_Betas(electrode,:,:,:));
                    for i=1:size(H1,3); E(:,i) = C*squeeze(H1(:,:,i))';end
                    lo = round(LIMO.design.nboot.*p./2);
                    hi = LIMO.design.nboot - lo;
                    sort_boot_betas = sort(E,2); % sort boot
                    b = sort_boot_betas(:,lo+1);
                    c = sort_boot_betas(:,hi);
                    if electrode == 1
                        if length(LIMO.design.electrode) == 1
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI optimized electrode', (1-p).*100);
                        end
                    else
                        mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                    end
                    
                else
                    
                    load Yr; load Res;
                    df = sum((Yr ~= NaN),3)-LIMO.design.nb_conditions; t = tcdf(1-p,df);
                    sigma2 = nansum((Res.^2 ./ repmat(df,[1 1 size(Res,3)])),3); clear Yr Res
                    X = LIMO.design.X(:,1:end-1); v = t.*sqrt(sigma2 ./ norm(X).^2);
                    b = Effect - v(electrode,:);
                    c = Effect + v(electrode,:);
                    if electrode == 1
                        if length(LIMO.design.electrode) == 1
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% CI optimized electrode', (1-p).*100);
                        end
                    else
                        mytitle = sprintf('Mean parameter difference between groups \n and %g%% CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                    end
                end
                
                plot(timevect,Effect,'LineWidth',3);
                hold on; plot(timevect,b,'--r','LineWidth',3);
                plot(timevect,c,'--r','LineWidth',3);
                xlabel('Time in ms','FontSize',16)
                ylabel('Amplitude (arbitrary unit)','FontSize',16)
                grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                hold on; plot(timevect,sig.*v(3),'r.','LineWidth',3);
                
                plotted_data = NaN(4,size(Effect,2));
                plotted_data(1,:) = c;
                plotted_data(2,:) = Effect;
                plotted_data(3,:) = b;
                plotted_data(4,:) = sig;
                assignin('base','Plotted_data',squeeze(plotted_data));
                
            else % plot all parameters
                
                if robust ==1
                    
                    index = 1;
                    [M, mask, mytitle] = limo_stat_values(Type,'One_way_variance_analysis',p,MCC,LIMO,choice,[]);
                    data = NaN(size(Betas,1),size(Betas,2),3*LIMO.design.nb_conditions);
                    
                    for effect_nb = 1:LIMO.design.nb_conditions
                        plotting_index(effect_nb) = index;
                        data(:,:,index+1) = Betas(:,:,effect_nb);
                        boot_data = squeeze(boot_H1_Betas(:,:,effect_nb,:));
                        if numel(size(boot_data)) == 2  % if one electrode only
                            boot_data = ones(1,size(boot_data,1),size(boot_data,2));
                            boot_data(1,:,:) = squeeze(boot_H1_Betas(:,:,effect_nb,:));
                        end
                        boot_data = sort(boot_data,3);
                        low = round((p*LIMO.design.nboot) / (2*size(LIMO.design.X,2)));
                        high = LIMO.design.nboot - low;
                        data(:,:,index) = boot_data(:,:,low);
                        data(:,:,index+2) = boot_data(:,:,high);
                        index = index+3;
                    end
                    
                    if electrode ==1 && length(LIMO.design.electrode) ~=1 % ie this is an optimized electrode
                        mytitle = {['Beta coefficient and ' num2str((1-p).*100) '% bootstrap CI'];['optimized electrode']};
                    else
                        mytitle = {['Beta coefficient and ' num2str((1-p).*100) '% bootstrap CI'];['at electrode ' LIMO.data.chanlocs(electrode).labels]};
                    end
                    
                else
                    
                    load Yr; load Res;
                    df = sum((Yr ~= NaN),3)-2; t = tcdf(1-p,df);
                    sigma2 = nansum((Res.^2 ./ repmat(df,[1 1 size(Res,3)])),3); clear Yr Res
                    
                    index = 1;
                    [M, mask, mytitle] = limo_stat_values(Type,'One_way_variance_analysis',p,MCC,LIMO,choice,[]);
                    data = NaN(size(Betas,1),size(Betas,2),3*LIMO.design.nb_conditions);
                    
                    for effect_nb = 1:LIMO.design.nb_conditions
                        plotting_index(effect_nb) = index;
                        X = LIMO.design.X(:,effect_nb);
                        data(:,:,index+1) = Betas(:,:,effect_nb);
                        data(:,:,index) = data(:,:,index+1) - t.*sqrt(sigma2 ./ norm(X).^2);
                        data(:,:,index+2) = data(:,:,index+1) + t.*sqrt(sigma2 ./ norm(X).^2);
                        index = index+3;
                    end
                    
                    if electrode ==1 && length(LIMO.design.electrode) ~=1 % ie this is an optimized electrode
                        mytitle = {['Beta coefficient and ' num2str((1-p).*100) '% CI'];['optimized electrode']};
                    else
                        mytitle = {['Beta coefficient and ' num2str((1-p).*100) '% CI'];['at electrode ' LIMO.data.chanlocs(electrode).labels]};
                    end
                end
                
                % figure
                figure;set(gcf,'Color','w'); hold on
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                plot(timevect,squeeze(data(electrode,:,plotting_index+1)),'LineWidth',3);
                plot(timevect,squeeze(data(electrode,:,plotting_index)),'--','LineWidth',2);
                plot(timevect,squeeze(data(electrode,:,plotting_index+2)),'--','LineWidth',2)
                xlabel('Time in ms','FontSize',16)
                ylabel('Amplitude (arbitrary unit)','FontSize',16)
                grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                plot(timevect,sig.*v(3),'r.','LineWidth',3);
                set(gca,'FontSize',14);
                
                plotted_data = NaN(size(data,2),size(data,3)+1);
                plotted_data(:,end) =sig; plotted_data(:,1:end-1) = squeeze(data(electrode,:,:));
                assignin('base','Plotted_data',squeeze(plotted_data));  % with the squeeze if N=1, reduces the matrix
                
            end
            
            
        elseif strcmp(LIMO.design.name,'One-way ANCOVA one electrode') || strcmp(LIMO.design.name,'One-way ANCOVA all electrodes')
            % --------------------------------------------------------------------------------------------------------------------------
            
            type = questdlg('2 types of plots are available','plotting option','plot main effect','plotting covariate(s)','plot main effect');
            load('Betas.mat');
            load('One_way_covariance_analysis')
            
            % check electrode
            [e,f,d]=size(Betas);
            if e > 1
                electrode = inputdlg('which electrode to plot','Plotting option');
                if isempty(electrode) || strcmp(cell2mat(electrode),'')
                    clear electrode; [v,e] = max(N_way_covariance_analysis(:,:,1)); [v,c]=max(v); electrode = e(c);
                else
                    electrode = eval(cell2mat(electrode));
                    if length(electrode) > 1
                        error('1 electrode only can be plotted')
                    elseif electrode > e
                        error('electrode number invalid')
                    end
                end
            else
                if length(LIMO.design.electrode) == 1
                    electrode = LIMO.design.electrode;
                else
                    electrode = 1;  % accomodates the fact that all matrices have the electrode dim (even = 1)
                end
            end
            
            
            % check bootstrap
            try
                load('boot_H1_Betas')
                if MCC ==1  %% if not corrected but boot data are there ask
                    choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
                    if strcmp(choice,'use empirical p values')
                        robust = 1;
                    else
                        robust = 0;
                    end
                else % MCC choosen and boot data are there
                    robust =1; choice = 'use empirical p values';
                end
            catch ME
                robust = 0;choice = 'use theoretical p values';
            end
            
            
            % plot
            if strcmp(type,'plot main effect')
                
                figure;set(gcf,'Color','w')
                [M, mask, mytitle] = limo_stat_values(Type,'One_way_covariance_analysis',p,MCC,LIMO,choice,1);
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                c = [eye(LIMO.design.nb_conditions-1) ones(LIMO.design.nb_conditions-1,1).*-1];
                C = limo_contrast_checking(LIMO.dir, LIMO.design.X, c);
                Effect = C*squeeze(Betas(electrode,:,:))';
                
                if robust == 1
                    H1 = squeeze(boot_H1_Betas(electrode,:,:,:));
                    for i=1:size(H1,3); E(:,i) = C*squeeze(H1(:,:,i))';end
                    lo = round(LIMO.design.nboot.*p./2);
                    hi = LIMO.design.nboot - lo;
                    sort_boot_betas = sort(E,2); % sort boot
                    b = sort_boot_betas(:,lo+1);
                    c = sort_boot_betas(:,hi);
                    if electrode == 1
                        if length(LIMO.design.electrode) == 1
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI optimized electrode', (1-p).*100);
                        end
                    else
                        mytitle = sprintf('Mean parameter difference between groups \n and %g%% bootstrap CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                    end
                    
                else
                    
                    load Yr; load Res;
                    df = sum((Yr ~= NaN),3)-LIMO.design.nb_conditions; t = tcdf(1-p,df);
                    sigma2 = nansum((Res.^2 ./ repmat(df,[1 1 size(Res,3)])),3); clear Yr Res
                    X = LIMO.design.X(:,1:LIMO.design.nb_conditions); v = t.*sqrt(sigma2 ./ norm(X).^2);
                    b = Effect - v(electrode,:);
                    c = Effect + v(electrode,:);
                    if electrode == 1
                        if length(LIMO.design.electrode) == 1
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Mean parameter difference between groups \n and %g%% CI optimized electrode', (1-p).*100);
                        end
                    else
                        mytitle = sprintf('Mean parameter difference between groups \n and %g%% CI at electrode %s (%g)', (1-p).*100, LIMO.data.chanlocs(electrode).labels,electrode);
                    end
                end
                
                plot(timevect,Effect,'LineWidth',3);
                hold on; plot(timevect,b,'--r','LineWidth',3);
                plot(timevect,c,'--r','LineWidth',3);
                xlabel('Time in ms','FontSize',16)
                ylabel('Amplitude (arbitrary unit)','FontSize',16)
                grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                sig = single(squeeze(mask(electrode,:))); sig(sig == 0) = NaN;
                hold on; plot(timevect,sig.*v(3),'r.','LineWidth',3);
                
                plotted_data = NaN(4,size(Effect,2));
                plotted_data(1,:) = c;
                plotted_data(2,:) = Effect;
                plotted_data(3,:) = b;
                plotted_data(4,:) = sig;
                assignin('base','Plotted_data',squeeze(plotted_data));
                
                
            else
                
                % check with regressor to plot
                if LIMO.design.nb_continuous
                    effect_nb = 3;
                else
                    limo_review(LIMO);
                    effect_nb = eval(cell2mat(inputdlg('which regressor to plot','Plotting option')));
                    if length(effect_nb)>1
                        go =0;
                        while go==0
                            errordlg('please select only one regressor');
                            effect_nb = eval(cell2mat(inputdlg('which covariate to plot','Plotting option')));
                            if length(effect_nb) == 1 ; go =1; end
                        end
                    end
                    close('Review Design')
                end
                
                % comnpute CI and plot
                data = NaN(size(Betas,1),size(Betas,2),3);
                if robust ==1
                    data(:,:,2) = Betas(:,:,effect_nb);
                    boot_data = squeeze(boot_H1_Betas(:,:,effect_nb,:));
                    if numel(size(boot_data)) == 2  % if one electrode only
                        boot_data = ones(1,size(boot_data,1),size(boot_data,2));
                        boot_data(1,:,:) = squeeze(boot_H1_Betas(:,:,effect_nb,:));
                    end
                    boot_data = sort(boot_data,3);
                    low = round((p*LIMO.design.nboot) / (2*size(LIMO.design.X,2)));
                    high = LIMO.design.nboot - low;
                    data(:,:,1) = boot_data(:,:,low);
                    data(:,:,3) = boot_data(:,:,high);
                    
                    % figure
                    figure;set(gcf,'Color','w')
                    timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                    plot(timevect,squeeze(data(electrode,:,2)),'LineWidth',3);
                    hold on; plot(timevect,squeeze(data(electrode,:,1)),'--r','LineWidth',2);
                    plot(timevect,squeeze(data(electrode,:,3)),'--r','LineWidth',2)
                    xlabel('Time in ms','FontSize',16)
                    ylabel('Amplitude (arbitrary unit)','FontSize',16)
                    if electrode ==1 && length(LIMO.design.electrode) ~=1 % ie this is an optimized electrode
                        mytitle = {['Beta coefficient ' num2str(effect_nb) ' and ' num2str((1-p).*100) '% bootstrap CI'];['optimized electrode']};
                    else
                        mytitle = {['Beta coefficient ' num2str(effect_nb) ' and ' num2str((1-p).*100) '% bootstrap CI'];['at electrode ' LIMO.data.chanlocs(electrode).labels]};
                    end
                    grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                    v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                    set(gca,'FontSize',14);
                    
                    plotted_data = NaN(size(data,2),3);
                    plotted_data(:,1) = squeeze(data(electrode,:,1)); % upper CI
                    plotted_data(:,2) = squeeze(data(electrode,:,2)); % data
                    plotted_data(:,3) = squeeze(data(electrode,:,3)); % lower CI
                    assignin('base','Plotted_data',squeeze(plotted_data));  % with the squeeze if N=1, reduces the matrix
                    
                else
                    
                    load Yr; load Res;
                    df = sum((Yr ~= NaN),3)-2; t = tcdf(1-p,df);
                    sigma2 = nansum((Res.^2 ./ repmat(df,[1 1 size(Res,3)])),3); clear Yr Res
                    if effect_nb > LIMO.design.nb_conditions
                        X = zscore(LIMO.design.X(:,effect_nb));
                    else
                        X = LIMO.design.X(:,effect_nb);
                    end
                    data(:,:,2) = Betas(:,:,effect_nb);
                    data(:,:,1) = data(:,:,2) - t.*sqrt(sigma2 ./ norm(X).^2);
                    data(:,:,3) = data(:,:,2) + t.*sqrt(sigma2 ./ norm(X).^2);
                    
                    % figure
                    figure;set(gcf,'Color','w')
                    timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                    plot(timevect,data(electrode,:,2),'LineWidth',3);
                    hold on; plot(timevect,data(electrode,:,1),'--r','LineWidth',2);
                    plot(timevect,data(electrode,:,3),'--r','LineWidth',2)
                    xlabel('Time in ms','FontSize',16)
                    ylabel('Amplitude (arbitrary unit)','FontSize',16)
                    if electrode ==1 && length(LIMO.design.electrode) ~=1 % ie this is an optimized electrode
                        mytitle = {['Beta coefficient ' num2str(effect_nb) ' and ' num2str((1-p).*100) '% bootstrap CI'];['optimized electrode']};
                    else
                        mytitle = {['Beta coefficient ' num2str(effect_nb) ' and ' num2str((1-p).*100) '% bootstrap CI'];['at electrode ' LIMO.data.chanlocs(electrode).labels]};
                    end
                    grid on; axis tight; title(mytitle,'FontSize',20); drawnow;
                    v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                    set(gca,'FontSize',14);
                    
                    plotted_data = NaN(size(data,2),3);
                    plotted_data(:,1) = squeeze(data(electrode,:,1)); % upper CI
                    plotted_data(:,2) = squeeze(data(electrode,:,2)); % data
                    plotted_data(:,3) = squeeze(data(electrode,:,3)); % lower CI
                    assignin('base','Plotted_data',squeeze(plotted_data));
                end
            end
            
        else
            errordlg('this file is not supported for this kind of plot','Nothing plotted')
        end
    end % closes type
end % closes if LIMO.level
end % closes the function

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                                   ROUTINES
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% color map
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function color_images_(scale,LIMO)

cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
set(gca,'XMinorTick','on','LineWidth',2)
try
    set(gca,'YTick',1:length(LIMO.data.expected_chanlocs));
catch ME
    set(gca,'YTick',1:length(LIMO.data.chanlocs));
end

ylabel('Electrodes','FontSize',14);
xlabel('Time in ms','FontSize',14);

if LIMO.Level == 1
    for i = 1 : length(LIMO.data.chanlocs)
        try
            label_electrodes{i} = LIMO.data.expected_chanlocs(i).labels;
        catch ME
            label_electrodes{i} = LIMO.data.chanlocs(i).labels;
        end
    end
    
else
    if isempty(LIMO.design.electrode)
        for i = 1 : length(LIMO.data.chanlocs)
            label_electrodes{i} = LIMO.data.chanlocs(i).labels;
        end
    else
        if length(LIMO.design.electrode) == 1
            label_electrodes = LIMO.design.electrode;
        else
            label_electrodes = ' ';
            ylabel('optimized electrode','FontSize',14);
        end
    end
end
set(gca,'YTickLabel', label_electrodes); 

%  OPTION TO ADJUST ELECTRODE LABEL SIZE --> DOES NOT WORK
%  IF SOMEONE KNOWS HOW TO DO IT PLEASE CONTACT US

%   for i = 1 : length(LIMO.data.chanlocs)
%       h = findobj('string',label_electrodes{i});
%       set(h,'FontSize',8)
%   end

% figure; axes;
% xlabel('not to change')
% ylabel('\bfElectrode channel')
% h = findobj('string','\bfElectrode channel');
% set(h,'FontSize',18)

end


%% time vector and label
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [timevect, label] = labels_time_(LIMO,ah)
interval = 50/(1000/LIMO.data.sampling_rate); % in frame
timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
zero_column = find(timevect == 0);
if isempty(zero_column) == 1 % in case it does not encompasses 0
    zero_column = 1;
end

if  LIMO.data.start < 0
    positive_label = timevect(zero_column:interval:end);
    negative_label = timevect(zero_column:-interval:1);
    if negative_label == 0
        negative_label = LIMO.data.start*1000;
    end
    label = [fliplr(negative_label(2:end)) positive_label];
else
    label = timevect(zero_column:interval:end);
end
set(ah,'XTick',find(timevect==label(1)):interval:find(timevect==label(end)),'XTickLabel',round(label));
end

