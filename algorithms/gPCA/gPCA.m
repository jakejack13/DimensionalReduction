function [mu,U,PCs] = gPCA(GPCAinput,ncol,timevec,k,M,max_iters,normalize,family)
%GPCA Takes matrix input and outputs reduced matrix using runGPCAv3.m 
%and rGPCAv3.r

basedir = '.'; % the directory where you store the runGPCAv3.R file
savedir = '.'; % where to write the mat file output

%save('tempGPCA.mat', 'GPCAinput', 'ncol', 'timevec', 'k', 'M', 'max_iters', 'normalize', 'family', '-v6')
%system(['Rscript ' basedir '/runGPCAv3.R ' savedir ' tempGPCA.mat GPCAoutput.mat']);
%system('rm tempGPCA.mat');

run runGPCAv3;

load('GPCAoutput.mat');
mu = GPCAoutput.mu;
U = GPCAoutput.U;
PCs = GPCAoutput.PCs;
system('rm GPCAoutput.mat');
end
