%% sample code to use generalized PCA with matlab
%
% first, install R and the relevant toolboxes:
%
% https://docs.google.com/document/d/1cD1KR9TdffHlH4N1Vjyn2e5cg0pStkydTQQGrslNpck/edit?usp=sharing
% 
% set up the companion file to this one, runGPCAv3.R
%
% after you run your gPCA, it will create a matfile called GPCAoutput.mat
%
% this contains the matrices PCs and U, as well as the inputs.
%
% The standard PCA I'm familiar with is X = U*S*V'
% 
% Where U and V are orthonormal, and S is a diagonal matrix containing the eigenvalues.
%
% In this case, assuming you use the log link function, it's X = exp(PCs*U')
%
% U is orthonormal, and PCs is only orthogonal. If you want to retrieve the eigenvalues, you can use
%
% S = sqrt(diag(PCs'*PCs))
%
% and you can get the orthonormal PCs by dividing each column of PCs by the corresponding entry in S.
%
% then you'll have X = exp(newPCs * diag(S) * U')
%
% Luke Sjulson, 2016-10-28







% you'll have to edit these
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
save('tempGPCA.mat', 'GPCAinput', 'ncol', 'timevec', 'k', 'M', 'max_iters', 'normalize', 'family', '-v6')

system(['Rscript ' basedir '/runGPCAv3.R ' savedir ' tempGPCA.mat GPCAoutput.mat']);

system('rm tempGPCA.mat');