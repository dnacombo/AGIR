function mire

Screen('Preference', 'VBLTimestampingMode', 1);
Screen('Preference','SkipSyncTests', 1);
Screen('Preference','VisualDebugLevel', 0);
Screen('Preference', 'SuppressAllWarnings', 1);
window = Screen(max(Screen('Screens')),'OpenWindow');
screen(window,'fillrect',128);
pathname = fileparts(which(mfilename));
img = imread([pathname filesep 'mire.tif']);
screen(window,'putimage',img);
Screen('Flip',window);

while 1
    if KbCheck
        while KbCheck; end
        break
    end
end
screen('Close',window)
