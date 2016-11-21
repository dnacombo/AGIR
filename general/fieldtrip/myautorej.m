%% myautorej
% based on http://bishoptechbits.blogspot.com/2011/05/automated-removal-of-independent.html

% note that this function is using some customised eeglab functions
% so we include these folders atop the path.
addpath('/home/chaumonm/Dropbox/MATLAB/general/eeglab');
addpath('/home/chaumonm/Dropbox/MATLAB/general');

if isempty(EEG.icawinv)
    error('No ica weights! Compute ICA on your data first.')
end
% here you select the methods you want to use to select your components.
rejects = [
    1 % autocorr
    1 % focal comp
    0 % SNR
    1 % trial variability
    1 % focal trial activity
    1 % Correlation with EOG
    ];

subplots = true; % if true, we'll make subplots of each of the measures. if false, we'll do one figure for each measure.
if subplots
    figure(24654654);
    set(gcf,'numbertitle', 'off','name','Auto component rejection measures')
end

NbMin = 1; % minimum number of the above rejection methods
% that must be satisfied to reject a component.

rejfields = {'icarejcorr' 'Autocorrelation' [         0         0    1.0000]
    'icarejfoc' 'Focal components' [         0    0.5000         0]
    'icarejSNR' 'Signal to noise ' [    0.8000         0         0]
    'icarejvar' 'Trial variability' [         0    0.7500    0.7500]
    'icarejfoctri' 'Focal trial variability' [    0.7500         0    0.7500]
    'icarejeogcorr' 'Correlation with EOG' [    0.7500    0.7500         0]
    };

ncomp= size(EEG.icawinv,2); % ncomp is number of components
icaacts = eeg_getdatact(EEG,'component',1:ncomp);
for ifield = 1:size(rejfields,1)
    EEG.reject.(rejfields{ifield}) = false(1,ncomp);
    EEG.reject.([rejfields{ifield} 'col']) = rejfields{ifield,3};
end

if rejects(1)
    %% Autocorrelation
    % Identifying noisy components
    %----------------------------------------------------------------
    dropautocorr=.5; % Will drop components with autocorrelation less than this value;
    mycorrint=round(12/(1000/EEG.srate)); %  find N pts corresponding to 12 ms
    rej = false(1,ncomp);
    for k=1:ncomp
        y=icaacts(k,:,:);
        yy=xcorr(mean(y,3),mycorrint,'coeff'); %autocorrel for pts sep by 12 ms
        autocorr(k) = yy(1);
        if yy(1) < dropautocorr
            rej(k)=true; % codes reject/accept for this component
        end
    end
    %----------------------------------------------------------------
    if subplots
        subplot(2,3,1);cla
    else
        figure(641984);clf
        set(gcf,'numbertitle','off','name','Autocorrelation')
    end
    set(gca,'fontsize',16)
    plot(autocorr);
    
    hold on
    xlim([0 ncomp+1]);
    s = std(autocorr);
    m = mean(autocorr);
    yl = ylim;xl = xlim;
    [x,y] = meshgrid(xl(1):.1:xl(2),yl(1):.1:yl(2));
    galpha = 1./(s*(2*pi)^.5).*exp(-(y-m).^2./(2.*s^2));
%     h = surf(x,y,-ones(size(y)));shading flat
%     color = [ 0 0 0]';
%     C = repmat(color,[1,size(y)]);
%     C = permute(C,[2 3 1]);
%     set(h,'alphadata',1-galpha,'alphadatamapping','scaled','facealpha','interp',...
%         'CData',C,'CDataMapping','direct')
    plot(xl,[dropautocorr dropautocorr],'r');
    
    plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{1,3},'markersize',40)
    title(['Rejected components, based on autocorrelation at 12 ms.'])
    autocorr(autocorr > dropautocorr) = NaN;
    plot(autocorr,'or')
    datacursormode on
    
    EEG.reject.icarejcorr = logical(rej);
    
end
if rejects(2)
    %% Focal activity
    %----------------------------------------------------------------
    focalICAout = 7; % zscore cutoff
    rej = false(1,ncomp);
    clear mywt
    for k=1:ncomp
        mywt(:,k) = sort(abs(zscore(EEG.icawinv(:,k))),'descend'); %sorts standardized weights in descending order
        if mywt(1,k) > focalICAout
            rej(k)=true;
        end
    end
    %----------------------------------------------------------------
    if subplots
        subplot(2,3,2);cla
    else
        figure(35165);clf
        set(gcf,'numbertitle','off','name','Focal activity')
    end
    set(gca,'fontsize',16)
    surf(mywt');
    xlim([1 ncomp+1])
    xl = xlim;
    ylim([1 ncomp+1])
    yl = ylim;
    zl = zlim;
    view(68,24);
    shading flat
    toplot = mywt(1,:);
    toplot(toplot < focalICAout) = NaN;
    hold on
    plot3([0 0],ylim,[focalICAout focalICAout]);
    plot3(ones(1,numel(toplot)),1:ncomp,toplot,'or');

    plot3(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,zl(2)-diff(zl)/20,'marker','.','color',rejfields{2,3},'markersize',40)

    xlabel('Sorted channels')
    ylabel('Components')
    zlabel('Standardized weights')
    title('Components with focal activity')
    
    EEG.reject.icarejfoc = logical(rej);
    datacursormode on
    
end

if rejects(3)
    %% Low Signal to noise components
    POI = [0 600];
    BL = [-1000 0];
    rejfields{3,2} = ['Signal to noise Time of interest ' num2str(POI,'%g ') ' and Baseline ' num2str(BL,'%g ') ' ms.'];
    %----------------------------------------------------------------
    %NB need to specify epochpt1 and epochpt2 to focus on period of interest (POI)
    SNRcut = 1; % default value for cutoff on ratio of SD post/pre baseline
    
    POIpts = timepts(POI);
    BLpts = timepts(BL);
    
    zz = zscore(icaacts,[],2);% zscore along time
    av1 = mean(zz(:,POIpts,:),3); % average activity in POI acros trials
    av2 = mean(zz(:,BLpts,:),3); % activity in baseline acros trials
    SNR = std(av1,[],2)./std(av2,[],2); % ratio of the standard deviations of activity and baseline
    rej = SNR < SNRcut;

    %----------------------------------------------------------------
    if subplots
        subplot(2,3,3);cla
    else
        figure(872);clf
        set(gcf,'numbertitle','off','name','Signal to Noise Ratio')
    end
    set(gca,'fontsize',16)
    plot(SNR);
    hold on
    xlim([0 ncomp+1]);
    xl = xlim; yl = ylim;
    plot(xl,[SNRcut SNRcut],'r');
    title(['Rejected components, based on signal to noise ratio between \newline Time of interest ' num2str(POI,'%g ') ' and Baseline ' num2str(BL,'%g ') ' ms.'])
    SNR(SNR > SNRcut) = NaN;
    plot(SNR,'or')
    plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{3,3},'markersize',40)
    datacursormode on
    
    EEG.reject.icarejSNR = rej';
    %----------------------------------------------------------------
end

if rejects(4)
    %% Trial variability (stdev)
    %----------------------------------------------------------------
    POI = [0 600];
    BL = [-500 -100];
    
    POIpts = timepts(POI);
    BLpts = timepts(BL);
    zz = zscore(icaacts,[],3);% zscore across trials
    B = zz(:,POIpts,:);
    meanabs = mean(abs(B),2); %average of the absolute value along time.
    stdlist = std(meanabs,[],3); %will be high if there is high trial-to-trial variation
    m = mean(stdlist);
    st = std(stdlist);
    
    crit=m+2*st; % criterion
    rej = stdlist' > crit;
    %----------------------------------------------------------------
    if subplots
        subplot(2,3,4);cla
    else
        figure(2154);clf
        set(gcf,'numbertitle','off','name','Trial variability (zscore)')
    end
    set(gca,'fontsize',16)
    plot(stdlist);
    hold on
    xlim([0 ncomp+1]);
    xl = xlim;yl = ylim;
    plot(xl,[crit crit],'r');
    title(['Rejected components, based on intertrial variability'])
    stdlist(stdlist < crit) = NaN;
    plot(stdlist,'or')
    plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{4,3},'markersize',40)
    datacursormode on
    
    EEG.reject.icarejvar = rej;
    %----------------------------------------------------------------
end
if rejects(5)
    %% Focal trial activity
    % Find components with focal trial activity (those that have activity
    % on just a few trials and are almost zero on others)
    %----------------------------------------------------------------

    focalICAout = 7; % zscore cutoff
    myact =sort(abs(zscore(mean(icaacts,2),[],3)),3,'descend'); % sorts standardized mean trial activity 
    % in descending order
    rej = myact(:,:,1) > focalICAout;
    EEG.reject.icarejfoctri = rej';
    
    %----------------------------------------------------------------
    if subplots
        subplot(2,3,5);cla
    else
        figure(88);clf
        set(gcf,'numbertitle','off','name','Focal trial activity')
    end
    set(gca,'fontsize',16)
    surf(squeeze(double(myact)));
    view(68,24);
    shading flat
    hold on
    ylim([0 ncomp+1])
    plot3([0 0],ylim,[focalICAout focalICAout]);
    
    toplot = myact(:,:,1);
    toplot(toplot < focalICAout) = NaN;
    plot3(ones(1,size(toplot)),1:ncomp,toplot,'or');
    xl = xlim; yl = ylim;zl = zlim;
    xlabel('Sorted trial activity')
    ylabel('Components')
    zlabel('Standardized mean trial activity')
    plot3(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,zl(2)-diff(zl)/20,'marker','.','color',rejfields{5,3},'markersize',40)
    datacursormode on
    
    title(['Rejected components, based focal trial activity'])
    %----------------------------------------------------------------
end
if rejects(6)
    %% Correlation with EOG
    corthreshV = .4;
    corthreshH = .3;
    
    Veogchannames = {'EXG1'};
    Heogchannames = {'EXG[23]'};
    Veogchan = chnb(Veogchannames);
    Heogchan = chnb(Heogchannames);
    VEOG = EEG.data(Veogchan,:,:);
    HEOG = EEG.data(Heogchan(1),:,:) - EEG.data(Heogchan(2),:,:);
    VEOG = VEOG(:);
    HEOG = HEOG(:);
    s = size(icaacts);
    ICs = reshape(icaacts,[ncomp prod(s(2:3))])';
    cV  = abs(corr(ICs,VEOG))';
    cH  = abs(corr(ICs,HEOG))';
    rejV = cV > corthreshV ;
    rejH = cH > corthreshH;
    
    %----------------------------------------------------------------
    if subplots
        subplot(2,3,6);cla
    else
        figure(21588);clf
        set(gcf,'numbertitle','off','name','Correlation with EOG')
    end
    set(gca,'fontsize',16)
    plot([cV;cH]');
    hold on
    xlim([0 ncomp+1]);
    xl = xlim;yl = ylim;
    plot(xlim,[corthreshV corthreshV ],'r:');
    plot(xlim,[corthreshH corthreshH ],'m:');
    
    title(['Rejected components, based on correlation with EOG'])
    ylabel('Correlation coef with EOG');
    xlabel('Component');
    cV(cV < corthreshV) = NaN;
    plot(cV,'or')
    cH(cH < corthreshH) = NaN;
    plot(cH,'om')
    plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{6,3},'markersize',40)
    datacursormode on
    
    EEG.reject.icarejeogcorr = [rejV|rejH];
    %----------------------------------------------------------------
end
  
%% Final computations
% combine in gcompreject field and pass to pop_selectcomps
EEG.reject.gcompreject = false(1,ncomp);
for ifield = 1:size(rejfields,1)
    EEG.reject.gcompreject = [EEG.reject.gcompreject ; EEG.reject.(rejfields{ifield})];
end
EEG.reject.gcompreject = sum(EEG.reject.gcompreject) >= NbMin;
% plotting
PLOTPERFIG = 35;
try
    close(hfig)
end
clear hfig
for ifig = 1:ceil((ncomp+1)/PLOTPERFIG)
    try
        pop_selectcomps(EEG, [1+(ifig-1)*PLOTPERFIG:min([ncomp+1,ifig*PLOTPERFIG])]);
    catch
        set(gca,'visible','off');
    end
    hfig(ifig) = gcf;
end;
%
ax = get(hfig,'children');
ax = vertcat(ax{:});
ax = ax(regexpcell(get(ax,'tag'),'topocomp'));
[dum idx] = sort(get(ax,'tag'));
ax = ax(idx);

% create markers next to each topoplot showing which threshold has been
% passed.
if not(numel(ax) == ncomp+1)
    error('can''t find all the topos for your components')
end
for i_comp = 1:ncomp
    axes(ax(i_comp))
    hold on
    for irej = 1:numel(rejects)
        if EEG.reject.(rejfields{irej,1})(i_comp)
            x = -.5 + (irej > 6);
            y = .5 - .1*irej-.3*(rem(irej-1,6)+1>3);
            scatter(x,y,'markerfacecolor',EEG.reject.([rejfields{irej} 'col']),'markeredgecolor',EEG.reject.([rejfields{irej} 'col']));
        end
    end
end

axes(ax(ncomp+1));
hold on
for irej = 1:numel(rejects)
    if rejects(irej)
        x = 0;
        y = .5 - .1*irej;
        
        scatter(x,y,'markerfacecolor',EEG.reject.([rejfields{irej} 'col']),'markeredgecolor',EEG.reject.([rejfields{irej} 'col']));
        text(x+.1,y,[rejfields{irej,2} ' (' num2str(sum(EEG.reject.(rejfields{irej,1}))) ')']);
    end
end
figure(hfig(1));
