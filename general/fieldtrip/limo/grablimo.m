function grablimo(GRAB)

% todo:
% 1: generate design matrix
% 2: align data for processing
% 3: run GLM
% 4: run 2dlevel...

%% so for each subject (1st dim of GRAB)
% we have n conditions along the other dimensions.
        LIMOCat = NaN * zeros(1,numel(EEG.event));
        
        LIMOCat(1,[YesNoise & Hits]) = 1;
        LIMOCat(1,[YesNoise & Misses]) = 2;
        LIMOCat(1,[NoNoise & Hits]) = 3;
        LIMOCat(1,[NoNoise & Misses]) = 4;
        LIMOCat(1,[YesNoise & FA]) = 5;
        LIMOCat(1,[YesNoise & CR]) = 6;
        LIMOCat(1,[NoNoise & FA]) = 7;
        LIMOCat(1,[NoNoise & CR]) = 8;
        save LIMOCat LIMOCat
