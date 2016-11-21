function mylimo_plots(expected_chanlocs,GRAB)

% allows multiple plots of parameter estimates
% most plots relies on median and quantile for
% robust estimates
%
% Cyril Pernet v1 25-May-2010
% -----------------------------
% Copyright (C) LIMO Team 2010

current_dir = pwd;
if nargin == 0
    error('expected chanlocs expected as input')
else
    if not(isstruct(expected_chanlocs))
        load(expected_chanlocs)
    end
end


for i = 1:numel(GRAB)
    [p f e] = fileparts(GRAB(i).name);
    Names{i} = [f e];
    Paths{i} = p;
    Files{i} = GRAB(i).name;
end


% check it's Betas.mat files and which param to test
% --------------------------------------------------
is_Betas = [];
for i=1:size(Names,2)
    if strcmp(Names{i},'Betas.mat')
        is_Betas(i) = 1;
    end
end

% if (isempty(is_Betas)) == 0 && sum(is_Betas) == size(Names,2)
%     parameters = eval(cell2mat(inputdlg('which parameters to test e.g [1:3]','parameters option')));
%     if isempty(parameters)
%         return
%     end
% else
%     errordlg('file selection failed, only Betas.mat files are supported'); return
% end


% match frames
% -------------
disp('matching frames across subjects')
first_frame = NaN(1,size(Paths,2));last_frame = first_frame;
start = last_frame; stop = start; sampling_rate = stop;
for i=1:size(Paths,2)
    cd (Paths{i});
    load LIMO;
    sampling_rate(i)          = LIMO.data.sampling_rate;
    first_frame(i)            = LIMO.data.trim1;
    last_frame(i)             = LIMO.data.trim2;
    start(i)                  = LIMO.data.start;
    stop(i)                   = LIMO.data.end;
    subj_chanlocs(i).chanlocs = LIMO.data.chanlocs;
    clear LIMO
end

if (sum(sampling_rate == sampling_rate(1))) ~= length(sampling_rate)
    errordlg('data have different sampling rates'); return
end

% put data together
%------------------

disp('gathering data ...'); index = 1;
for i=1:size(Paths,2) % for each subject
    fprintf('processing subject %g',i); disp(' ')
    cd(Paths{i});
    load(Names{i});
    begins_at = max(first_frame) - first_frame(i) + 1;
    ends_at = size(Betas,2) - (last_frame(i) - min(last_frame));
    
    for j=1:size(Betas,3);%parameters
        tmp =  squeeze(Betas(:,:,j));
        data(:,:,j,i) = limo_match_elec(subj_chanlocs(i).chanlocs,expected_chanlocs,begins_at,ends_at,tmp);
        clear tmp
    end
end


% Make various plots
% -------------------
cd(current_dir); go = 1;

while go == 1
    %try
    K = local_menu;
    
    switch K
        
        % correlation of parameters
        % ------------------------
        case 1
            
            p = eval(cell2mat(inputdlg('enter 1 parameter of interest','parameter selection')));
            if p > size(data,3)
                error('parameter value > data size');
            elseif length(p) > 1
                error('only one value allowed')
            end
            
            figure('Name',['Plot of median parameter ',num2str(p)]); set(gcf,'Color','w');
            surf(nanmedian(data(:,:,p,:),4)); axis tight, grid on; ylabel('electrodes') ; xlabel('frames');shading interp
            K = local_menu;
            
            % correlation of parameters
            % ------------------------
        case 2
            
            alpha = 5/100;
            if size(data,3) ~= 2
                p = inputdlg('enter 2 parameters for the correlation e.g. [1 2]','parameter selection');
                p = eval(['[' p{1} ']']);
            else
                p = [1 2];
            end
            
            if p > size(data,3)
                error('parameter value > data size');
            elseif length(p) > 2
                error('only two values allowed')
            end
            x = nanmedian(data(:,:,p(1),:),4);
            y = nanmedian(data(:,:,p(2),:),4);
            
            % for each frame across electrodes
            [r,pval] = corr(x,y); % get r and p and mask of significant values after bootsrap
            [pID,pN] = limo_FDR(pval(pval<0.5),.05); % correct for multiple comparisons
            figure('Name',['Correlation matrices of median parameters ',num2str(p)]);set(gcf,'Color','w');
            subplot(1,2,1); imagesc(r.*single(pval<pN)); % image r < p
            axis square; title({['correlations across electrodes']; ['with FDR correction']},'Fontsize',15);
            xlabel('time');ylabel('time')
            
            % for each electrode across frames
            [r,pval]=corr(x',y','type','Spearman'); % get r and p and mask of significant values after bootsrap
            [pID,pN] = limo_FDR(pval(pval<0.5),.05); % correct for multiple comparisons
            subplot(1,2,2); imagesc(r.*single(pval<pN)); % image r < p
            axis square; title({['correlations across time']; ['with FDR correction']},'Fontsize',15);
            xlabel('electrodes');ylabel('electodes')
            K = local_menu;
            
            % joint distribution and plots in time
            % -----------------------------------
        case 3
            
            % select parameters
            if size(data,3) ~= 2
                p = inputdlg('enter 2 parameters for joint histograms e.g. [1 2]','parameter selection');
                p = eval(['[' p{1} ']']);
            else
                p = [1 2];
            end
            
            if p > size(data,3)
                error('parameter value > data size');
            elseif length(p) > 2
                error('only two values allowed')
            end
            
            % select time window
            time_vector = round(max(start)*1000:(1000/sampling_rate(1)):min(stop)*1000);
            S = min(time_vector); E = max(time_vector);
            t = eval(cell2mat(inputdlg([' enter time vector e.g. [1:20:100] min ',num2str(S),' max ',num2str(E)],'timing information')));
            if sum(t<S)~=0 || sum(t>E)~=0
                error('wrong timing information');
            end
            
            for i=1:length(t)
                difference = time_vector - t(i);
                [v,position]=find(difference==0);
                if v == 1
                    frames(i)=position;
                else
                    difference = difference .* (difference > 0);
                    difference(difference == 0) = NaN;
                    [v,frames(i)]=min(difference);
                end
            end
            x = squeeze(data(:,frames,p(1),:));
            y = squeeze(data(:,frames,p(2),:));
            
            % make figure
            % ideally the axes would be labelled with time_info =time_vector(frames);
            figure('Name',['Joint scatter of median parameters ',num2str(p)]); set(gcf,'Color','w');
            subplot(2,2,1); plotmatrix(nanmedian(x,3)); xlabel('time'); ylabel('time'); title(sprintf('scatter plots per frame parameter %g',p(1)),'Fontsize',12);
            subplot(2,2,2); plotmatrix(nanmedian(x,3),nanmedian(y,3)); title('joint scatter plots of parameters','Fontsize',12);
            subplot(2,2,4); plotmatrix(nanmedian(y,3)); xlabel('time'); ylabel('time'); title(sprintf('scatter plots per frame parameter %g',p(2)),'Fontsize',12);
            K = local_menu;
            
            % classic box-plot
            % ----------------
        case 4
            
            [e,f,r,n]=size(data);
            if e>1
                electrode = eval(cell2mat(inputdlg('which electrode to plot','Plotting option')));
                if electrode > e
                    error('electrode number invalid')
                end
            end
            
            % select time window
            time_vector = round(max(start)*1000:(1000/sampling_rate(1)):min(stop)*1000);
            S = min(time_vector); E = max(time_vector);
            t = eval(cell2mat(inputdlg([' enter time vector e.g. [1:20:100] min ',num2str(S),' max ',num2str(E)],'timing information')));
            if sum(t<S)~=0 || sum(t>E)~=0
                error('wrong timing information');
            end
            
            for i=1:length(t)
                difference = time_vector - t(i);
                [v,position]=find(difference==0);
                if v == 1
                    frames(i)=position;
                else
                    difference = difference .* (difference > 0);
                    difference(difference == 0) = NaN;
                    [v,frames(i)]=min(difference);
                end
            end
            
            
            % select parameters
            if size(data,3) ~= 2
                p = inputdlg('enter 2 parameters for joint histograms e.g. [1 2]','parameter selection');
                p = eval(['[' p{1} ']']);
            else
                p = [1 2];
            end
            
            if p > size(data,3)
                error('parameter value > data size');
            elseif length(p) > 2
                error('only two values allowed')
            end
            
            % do the figures
            time_info =time_vector(frames);
            for i = 1:length(frames)
                figure('Name','Parameter boxplot');set(gcf,'Color','w')
                boxplot(squeeze(data(electrode,frames(i),p,:))','notch','on');
                xlabel('Regressor(s)','FontSize',16)
                mytitle = sprintf('Electrode %g time %g ms', electrode, time_info(i));
                grid on; title(mytitle,'FontSize',20); drawnow;
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                ylabel('parameter value (abstract unit)','FontSize',16)
                drawnow
            end
            
            K = local_menu;
        otherwise
            go = 0;
            break
    end
end
% end


%% menu
function K = local_menu

K = menu('Choose a type of plot','Surf ERP space','Correlation matrices','Joint scatter plots','Box Plots','Quit') ;


