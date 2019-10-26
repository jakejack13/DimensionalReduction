basedir = '.'; % the directory where you store the runGPCAv3.R file
savedir = '.'; % where to write the mat file output


%GPCAinput = round(5*rand(100, 8)); % generating a fake input

%timevec = 0:100;                    % just illustrating that you can include other info


% parameters to pass to R in the mat file
%ncol = size(GPCAinput, 2);
%k = ncol; % how many PCs to compute
%M = 4; % value to approximate the saturated model - default 4, can set to 2 if it doesn't converge
%max_iters = 1000; % maximum number of iterations - default 1000
%normalize = 0;  % whether or not to normalize so each neuron has the same weight. 0 for no, 1 for yes.
%family = 'poisson';  


% whatever variables you save to the temporary mat file get passed to R, which save them in its mat file along with the output
save('tempLPCA.mat', 'LPCAinput', 'ncol', 'timevec', 'k', 'M', 'max_iters', 'normalize', 'family', '-v6')

system(['Rscript ' basedir '/runLPCA.R ' savedir ' tempLPCA.mat LPCAoutput.mat']);

system('rm tempLPCA.mat');