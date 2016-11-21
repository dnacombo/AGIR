function list_trig = trig_theoretical(C)

% function giving the list of triggers that should be recorded during all
% the task.

%definitions
trigger_1 = 200; %% each heart beat
trigger_2 = 50; %% receing HB in Sync conditions
trigger_3 = 100; %% receing HB in ASync conditions
trigger_4 = 30; %% Asynchrone definition
trigger_5 = 68; %% Synchrone definition and begin scale
trigger_6 = 80; %% begining of each session
trigger_7 = 240; %% first trial for each block
trigger_8 = 112; %% flip passive
trigger_9 = 130; %% moteur run passive
trigger_10 = 150; %% flip active
trigger_11 = 165; %% motor run active
trigger_12 = 180; %% begin scale
trigger_13 = 250; %% end scale
trigger_14 = 220; %% end task
trigger_15 = 230; %% after last trial break while loop
trigger_16 = 10; %% question 1 : no
trigger_17 = 20; %% question 1 : yes
trigger_18 = 40; %% question 2 answer: respiration
trigger_19 = 60; %% question 2 answer: erreur
trigger_20 = 70; %% question 2 answer: r�activit�
trigger_21 = 90; %% question 2 answer: force
trigger_22 = 110; %% question 2 answer: coeur
trigger_23 = 120; %% question 2 answer: autre
trigger_24 = 140; %% start complementary test
trigger_25 = 160; %% answer interoception score, answer = Sync
trigger_26 = 170; %% answer interoception score, answer = ASync
trigger_27 = 190; %% end of training, begining interoception test.


% columns in the matrix of trials we created, to keep in mind
col_nb_total_trial=1;
col_sync=2;
col_percentage = 3;
col_this_trial_type = 4;
col_what_proportion_tells_type = 5;


%number of blocks in each session
number_of_blocks_per_sequence = 9;

list_trig = [];

for this_block = 1 : size(C, 2)
            total_trial = size(C{this_block}, col_nb_total_trial);

    
    if this_block == 1
        list_trig = [list_trig, trigger_6];
        
        if C{this_block}(1, col_sync) ==1;
            list_trig = [list_trig, trigger_5];
            
            
        elseif C{this_block}(1, col_sync) ==2;
            list_trig = [list_trig, trigger_4];
        end
    end
        
    if mod(this_block, number_of_blocks_per_sequence) ==  0; % for the stop session message
        list_trig = [list_trig, trigger_6];
    end    
    
    for this_trial = 1 : 10*(C{this_block}(col_nb_total_trial))
        
        if this_trial == 1
            list_trig = [list_trig, trigger_7];            
        end
                
        Trial_Type = C{this_block}(this_trial, col_this_trial_type);
        
        if C{this_block}(this_trial, col_what_proportion_tells_type) ==1;
            if Trial_Type == 1
                type_of_this_trial = 1; %active trial
            elseif Trial_Type == 2
                type_of_this_trial = 2; % passive trial
            end
        elseif  C{this_block}(this_trial, col_what_proportion_tells_type) ==2;
            if Trial_Type == 1
                type_of_this_trial = 2; % passive
            elseif Trial_Type == 2
                type_of_this_trial = 1; % active
            end
        end
    
    if type_of_this_trial == 1  %active
        list_trig = [list_trig, trigger_10];
        
    end
    
    if type_of_this_trial == 2  % passive
        list_trig = [list_trig, trigger_8];
        list_trig = [list_trig, trigger_9];
    end
    
    if this_trial == total_trial
        list_trig = [list_trig, trigger_15];        
    end
end
    
    
    %% scale
    if this_block < size(C, 2)
        
        if C{this_block + 1}(1, col_sync) ==1;
            list_trig = [list_trig, trigger_5];
        elseif C{this_block + 1}(1, col_sync) ==2;
            list_trig = [list_trig, trigger_4];
        end
    end
    
    list_trig = [list_trig, trigger_12];
    list_trig = [list_trig, trigger_13];
    
end