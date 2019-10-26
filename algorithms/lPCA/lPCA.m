function [mu,U,PCs] = lPCA(LPCAinput,ncol,timevec,k,M,max_iters,normalize,family)
%LPCA Takes matrix input and outputs reduced matrix using runLPCAv3.m 
%and runLPCAv3.R

basedir = '.'; % the directory where you store the runGPCAv3.R file
savedir = '.'; % where to write the mat file output

%save('tempGPCA.mat', 'GPCAinput', 'ncol', 'timevec', 'k', 'M', 'max_iters', 'normalize', 'family', '-v6')
%system(['Rscript ' basedir '/runGPCAv3.R ' savedir ' tempGPCA.mat GPCAoutput.mat']);
%system('rm tempGPCA.mat');

run runLPCA;

load('LPCAoutput.mat');
mu = LPCAoutput.mu;
U = LPCAoutput.U;
PCs = LPCAoutput.PCs;
system('rm LPCAoutput.mat');
end
