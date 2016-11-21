goto home

[dum myrev] = system('svn info eeglab13');
[dum currrev] = system('svn info https://sccn.ucsd.edu/svn/software/eeglab');
% myrev = regexp(myrev,'Revision: (\d*)','tokens');
% currrev = regexp(currrev,'Revision: (\d*)','tokens');
% myrev = myrev{1}{1};
% currrev = currrev{1}{1};


rep = questdlg({['Your version is '] myrev [' most recent version is '] [currrev] ' do you want to continue?'});

switch rep
    case 'Yes'
    system('svn checkout https://sccn.ucsd.edu/svn/software/eeglab eeglab13');
end        

goback