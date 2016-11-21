function [keys] = GetResponseKeys(reps,cmds,varargin)
    
if not(exist('cmds','var'))
    cmds = '';
end
if not(isempty(regexp(cmds,'kBoard', 'once' )))
    kBoard = varargin{1};
else
    kBoard = -1;
end

while KbCheck(kBoard);end
commandwindow;
keys = zeros(1,numel(reps));
for i_rep = 1:numel(reps)
    fprintf(['Press the ' reps{i_rep} ' key...']);
    while 1
        [ keyIsDown, seconds, keyCode ] = KbCheck(kBoard);
        if keyIsDown
            while KbCheck(kBoard); end
            keys(i_rep) = find(keyCode);
            fprintf([' (KeyCode = ' num2str(keys(i_rep)) ')'])
            break
        end
    end
    fprintf('\n');
end

fprintf('Thanks!\n');

if not(isempty(regexp(cmds,'struct', 'once' )))
    out = struct;
    for i = 1:numel(reps)
        out.(regexprep(reps{i},'[^a-zA-Z0-9]','')) = keys(i);
    end
    keys = out;
end    