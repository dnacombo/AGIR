function limo_1stLevel_addcontrasts(Contrasts,opts)

if not(exist('opts','var'))
    opts = struct;
end
def.rootdir = cd;
def.re = '.*.set';
def.rmchan = {'EXG.*'};
def.timewindow = [-Inf Inf];
opts = setdef(opts,def);

% we load all subjects' info.
cd(opts.rootdir)
GRAB = datagrabber(opts.re, 'loadmode','info');

for i_des = 1:numel(Contrasts)
    cd(opts.rootdir)
    for isubj = 1:numel(GRAB)
        LIMO = [];
        EEG = [];
        disp(fileparts(GRAB(isubj).name));
        
        sujdir = fileparts(GRAB(isubj).name);
        
        cd([sujdir filesep Contrasts(i_des).designname]);
        
        %%% load useful information and data.
        load LIMO
        load Yr
        load Betas
        
        f = GRAB(isubj).name;
        
        % we will remove contrast regressors that could be missing for
        % some subjects (due to zero trials...)
        missing_conds = find(diff([0; LIMO.design.conditions]) > 1);
        
        % clear any existing previous contrasts because we don't want to mixup
        % things.
        delete con_*.mat ess_*.mat
        LIMO.contrast = {};
        for i_c = 1:size(Contrasts(i_des).Cs,1)% then for each of them,
            % compute contrasts %%
            cs = Contrasts(i_des).Cs{i_c,1};
            missing_conds(missing_conds > size(cs,2)) = [];
            cs(:,missing_conds) = [];
            % the two following lines result in setting all non specified
            % regressors to zero in the contrast.
            ctrsts = zeros(size(cs,1),LIMO.design.nb_conditions + LIMO.design.nb_continuous + LIMO.design.intercept);
            ctrsts(:,1:size(cs,2)) = cs;
            
            % update LIMO structure
            LIMO.contrast{i_c}.C = ctrsts;
            LIMO.contrast{i_c}.V = Contrasts(i_des).Ctypes{i_c};
            %%% initialize saved files.
            % name the contrast
            switch LIMO.contrast{i_c}.V
                case 'T'
                    con = zeros(size(Yr,1),size(Yr,2),3); % 3 dims = C*Beta/t/p
                    filename = sprintf('con_%g.mat',i_c);
                    save ([filename], 'con'); clear con;
                case 'F'
                    ess = zeros(size(Yr,1),size(Yr,2),size(LIMO.contrast{i_c}.C,1)+2); % 3 dims = C*Beta/F/p
                    filename = sprintf('ess_%g.mat',i_c);
                    save ([filename], 'ess'); clear ess;
            end
            str = LIMO.contrast{i_c}.V;
            if size(Contrasts(i_des).Cs,2) == 1;
                for i = 1:size(Contrasts(i_des).Cs{i_c,1},2)
                    switch Contrasts(i_des).Cs{i_c,1}(i)
                        case 1
                            str = [str ' +(' LIMO.design.XNames{i} ')'];
                        case -1
                            str = [str ' -(' LIMO.design.XNames{i} ')'];
                        case 0
                        otherwise
                            str = [str ' +' num2str(Contrasts(i_des).Cs{i_c,1}(i),'%.2g') '(' LIMO.design.XNames{i} ')'];
                    end
                end
            else
                str = Contrasts(i_des).Cs{i_c,2};
            end
            LIMO.contrast{i_c}.Name = str;
            
            fprintf(['====== Contrast :\n'])
            disp(ctrsts);
            disp(LIMO.contrast{i_c}.Name);
            % starting computation
            fprintf('Analysing electrode ');
            str = '';
            for electrode = 1:size(Yr,1)
                % display electrode number
                for i = 1:numel(str)
                    fprintf('\b');
                end
                str = sprintf('%02g',electrode);
                fprintf(str);
                if electrode == size(Yr,1)
                    fprintf('\n');
                end
                %%% compute contrast.
                switch LIMO.contrast{i_c}.V
                    case 'T'
                        go = limo_contrast_checking(LIMO.contrast{i_c}.C,LIMO.design.X);
                        if not(go)
                            error('check your contrast')
                        end
                        con(electrode,:,:) = mylimo_contrast(squeeze(Yr(electrode,:,:))', squeeze(Betas(electrode,:,:))', LIMO, 'T',1,i_c);
                    case 'F'
                        go = limo_contrast_checking(LIMO.contrast{i_c}.C,LIMO.design.X);
                        if not(go)
                            error('check your contrast')
                        end
                        ess(electrode,:,:) = mylimo_contrast(squeeze(Yr(electrode,:,:))', squeeze(Betas(electrode,:,:))', LIMO, 'F',1,i_c);
                end
            end
            % and save.
            switch LIMO.contrast{i_c}.V
                case 'T'
                    filename = sprintf('con_%g.mat',i_c);
                    condesc = LIMO.contrast{i_c};
                    save ([filename], 'con','condesc'); clear con condesc
                case 'F'
                    filename = sprintf('ess_%g.mat',i_c);
                    condesc = LIMO.contrast{i_c};
                    save ([filename], 'ess','condesc'); clear ess
            end
        end
        save LIMO LIMO
        
        %%% saving in EEG format.
        EEG = eeg_emptyset;
        EEG.chanlocs = LIMO.data.chanlocs;
        EEG.srate = LIMO.data.sampling_rate;
        EEG.xmin = LIMO.data.start;
        EEGi = EEG;
        clear bigcon
        c = dir('con_*.mat');
        for i = 1:numel(c)
            load(c(i).name);
            EEG = EEGi;
            EEG.setname = ['Contrast (T) ' condesc.Name];
            EEG.comments = sprintf('Contrast: %s',fullfile(cd,c(i).name));
            
            EEG = mat2eeglab(con(:,:,2),EEG);% store the T value.
            bigcon(:,:,i) = con(:,:,2);
            pop_saveset(EEG, 'filename', strrep(c(i).name,'.mat','.set'),'filepath',cd);
        end
        clear con
        if numel(c) > 0
            EEG = EEGi;
            EEG.setname = ['Contrast (T) (all)'];
            EEG.comments = sprintf('TContrasts: all in different epochs');
            
            EEG = mat2eeglab(bigcon,EEG);% all the T values
            pop_saveset(EEG, 'filename', 'con.set','filepath',cd);
        end
        
        c = dir('ess_*.mat');
        for i = 1:numel(c)
            load(c(i).name);
            EEG = EEGi;
            EEG.setname = ['Contrast (F) ' condesc.Name];
            EEG.comments = sprintf('Contrast: %s',fullfile(cd,c(i).name));
            
            EEG = mat2eeglab(ess(:,:,end-1),EEG); % store the F value
            bigess(:,:,i) = ess(:,:,end-1);
            pop_saveset(EEG, 'filename', strrep(c(i).name,'.mat','.set'),'filepath',cd);
        end
        clear ess
        if numel(c) > 0
            EEG = EEGi;
            EEG.setname = ['Contrast (F) (all)'];
            EEG.comments = sprintf('FContrasts: all in different epochs');
            
            EEG = mat2eeglab(bigess,EEG);
            
            pop_saveset(EEG, 'filename', strrep(c(i).name,'.mat','.set'),'filepath',cd);
        end
    end
end





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
