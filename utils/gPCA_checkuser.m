function [baseDir] = gPCA_checkuser()
% this function returns the working base directory that is specific for 
% different computers (directory of repository)

[~, hostname] = system('hostname');

if contains(hostname, 'Luke-HP-laptop')
	baseDir = 'C:\Users\lukes\Desktop\gPCA';
elseif contains(hostname, 'freyja')
	baseDir = '/home/labadmin/Documents/MATLAB/VarTests/gPCA';
	
else 
	warning('your computer is not known - please add to gPCA_checkuser.m');
	baseDir = [];
end

