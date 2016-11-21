function limo_2dLevel(designs,opts)

if not(exist('opts','var'))
    opts = struct;
end
def.rootdir = cd;
def.expected_chanlocs = cdhome('/expected_chanlocs.mat');
def.chanlocs = cdhome('/EEGchanlocs.mat');
def.nboot = 200;
def.dirre = '(?<suj>S.*_strength)/';
opts = setdef(opts,def);

%% limo_2dlevel script
global LIMO EEG
LIMO = [];
EEG = [];
load(opts.expected_chanlocs)
load(opts.chanlocs)
ELECS = {chanlocs.labels};
c = regexpcell({expected_chanlocs.labels},ELECS,'exact');
expected_chanlocs = expected_chanlocs(c);
channeighbstructmat = channeighbstructmat(c,c);
clear c


% designs = {'Masto_Cont_ConsensusTargetPrime_Relatedness' 'Masto_Cat_ConsensusTarget_Relatedness' 'Masto_Cat_ExtremeConsensusTarget_Relatedness' 'Masto_Cat_ExtremeConsensusPrime_Relatedness'};
% {'Masto_Cat_ConsensusTarget_Relatedness' 'Masto_Cat_ConsensusPrime_Relatedness' 'Masto_Cont_Relatedness_Consensus' 'Masto_Cont_Orthog_Relatedness_Consensus'};
%
%{'Cat_HL_Relatedness'};
%%
for i_des = 1:numel(designs)
    cd(opts.rootdir)
    designname = designs{i_des};
    
    
    stuff2exclude = 'AllSuj|Res.set|Betas.set|Continuous.set|Condition_effect.set|con.set';
    re = [opts.dirre designname '/(?<param>.*).set'];
    GRAB = datagrabber(re,'exclude',stuff2exclude);%,'waitafterlisting',1);
    if not(isdir('AllSuj'))
        mkdir AllSuj
    end
    cd AllSuj
    if not(isdir(designname))
        mkdir(designname)
    end
    cd(designname);
    
    %%
    GB = GRAB;
    clear GRAB
    for i_p = 1:size(GB,2)
        GRAB = GB(:,i_p);
        disp(GRAB(1).param);
        [dum name] = fileparts(GRAB(1).name);
        humanname = GRAB(1).EEG.setname;
        y = single(NaN([numel(ELECS), size(GRAB(1).data,2), numel(unique([GRAB.sujidx]))]));
        for i_d = 1:numel(GRAB)
            if strcmp(GRAB(i_d).param,'Res')% if we're dealing with residuals, 
                GRAB(i_d).range.trials.oper = 'mean';
                GRAB(i_d) = datagrabber(GRAB(i_d),'range',GRAB(i_d).range,'rerange',true);
            end
            % this is about keeping consistent order of electrodes.
            yd = NaN(numel(ELECS),size(GRAB(i_d).data,2));
            idxELECS = regexpcell(ELECS,GRAB(i_d).range.channels.chans,'exact');
            idxelecs = regexpcell(GRAB(i_d).range.channels.chans,ELECS(idxELECS),'exact');
            yd(idxELECS,:) = GRAB(i_d).data(idxelecs,:);
            y(:,:,GRAB(i_d).sujidx) = yd;
        end
        gp = unique([GRAB.sujidx]');
        factor_levels = numel(unique([GRAB.cond]));
        % limo_random_robust(7,y,gp,factor_levels,opts.nboot);%limo_rep_anova(y,gp,factor_levels);
        [stats one_sample boot_one_sample] = mylimo_random_robust(1,y,1,opts.nboot);%limo_rep_anova(y,gp,factor_levels);
        save(['one_sample_' name],'one_sample');
        save(['boot_one_sample_' name],'boot_one_sample');
        
        save(name,'-struct','stats');%,'one_sample','boot_one_sample');
        clear boot_one_sample one_sample;
        %return
        %% copy .mat in an EEGset
        EEG = eeg_emptyset;
        EEG.chanlocs = chanlocs;
        EEG.srate = GRAB(1).EEG.srate;
        EEG.xmin = GRAB(1).EEG.xmin;
        EEG.xmax = GRAB(1).EEG.xmax;
        
        EEG.setname = ['2d level parameter ' humanname ];
        EEG.statsfile = [name '.mat'];
        
%         EEG = mat2eeglab(stats.t,EEG);
        EEG = mat2eeglab(stats.m,EEG);
        pop_saveset(EEG, 'filename', [name '.set'],'filepath',cd);
        
    end
end
disp done
