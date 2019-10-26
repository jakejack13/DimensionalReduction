function [] = pca_process_v1(dataDir,binF)

%PCA_PROCESS_V1 Processes NAcSpikes data and stores into mat file
%dataDir = directory where the NAcSpikes.mat file is stored

% Jacob Kerr, 2019

%% Parameters for functions

%Going from gPCA subdirectory to data directory to access data
[baseDir] = gPCA_checkuser;
cd(dataDir) 
load('NAcSpikes.mat')

%Setting up input settings
%U = [3,2,1; 2,-3 0; 0,0,0];
%S = [5,0,0;4,0,0;2,0,0];
%t = 0:0.01:4*pi;
%V = [sin(t) + 1; sin(1.7*t) + 1; sin(5.3*t) + 1];
%I = U*S*V;

%input = round(rand(100,6));

input = extract_data(NAcSpikes,binF); %binning data from raw

ncol = size(input,2);
timevec = 0:100;
k = ncol;
M = 4;
max_iters = 1000;
normalize = 0;
family = 'poisson';

%Going back to gPCA subdirectory to run the algorithms
cd(baseDir);
cd(['algorithms/gPCA']);


%Output legend
%X_PC = Output matrix (principle components)
%X_var = Variance of principle components (eigenvalues)
%X_vec = eigenvectors


%Algorithm execution

%PCA
tic;
[pca_PC,~,pca_var] = pca(input);
pca_time = toc;
%pca_recon = pca_vec * transpose(pca_PC);

%NMF
tic;
[nnmf_var,nnmf_PC] = nnmf(input, k);
nnmf_time = toc;
%nnmf_recon = nnmf_var * nnmf_PC;

%SVD
tic;
[~,svd_var,svd_PC] = svd(input,0);
svd_time = toc;

%GPCA
tic;
[~,gpca_PC,gpca_var] = gPCA(input, ncol, timevec, k, M, max_iters, normalize, family);
gpca_time = toc;
%gpca_recon = gpca_var * transpose(gpca_PC);

%EPCA
tic;
[epca_PC,epca_var,~,~,~,~,~,~,~,~] = exp_fam_pca(input,family,k,0,k);
epca_time = toc;
%epca_recon = epca_vec * transpose(epca_PC);
    
%Moving to lPCA directory
cd(['..']);
cd(['lPCA']);

%Rounding data to binomial distribution for lPCA
%WARNING: input data after this is rounded down and is not the original
%data
for idx = 1:numel(input)
    if input(idx) > 1
        input(idx) = 1;
    end
end

%LPCA
tic;
[~,lpca_PC,lpca_var] = lPCA(input, ncol, timevec, k, M, max_iters, normalize, family);
lpca_time = toc;
L = lpca_var * transpose(lpca_PC);
%lpca_recon = lpca_var * transpose(lpca_PC);

%Compilation of algo runtimes
algo_time = [pca_time;nnmf_time;gpca_time;lpca_time;epca_time];
%algo_time = [pca_time;nnmf_time;svd_time;gpca_time;lpca_time;epca_time];

%Moving back to base directory 
cd(baseDir);

%% Scaling tests

%Output legend
%X_scale = scaled eigenvalues (percentages of variance)
%X_var_norm = normalized eigenvalues


%Converting raw eigenvalues to percentages of variance
pca_var_scale = pca_var ./ sum(pca_var);

nnmf_var_norm = vecnorm(nnmf_var);
nnmf_var_scale = nnmf_var_norm ./ sum(nnmf_var_norm);

svd_var_scale = svd_var ./ sum(svd_var);

gpca_var_norm = vecnorm(gpca_var);
gpca_var_scale = gpca_var_norm ./ sum(gpca_var_norm);

lpca_var_norm = vecnorm(lpca_var);
lpca_var_scale = lpca_var_norm ./ sum(lpca_var_norm);

epca_var_scale = epca_var/sum(epca_var);

%Normalizing PC1 of PCA and scaling other matrices by this factor
var_scale = 1 ./ pca_var_scale(1);
pca_var_scale = transpose(pca_var_scale) .* var_scale;
nnmf_var_scale = nnmf_var_scale .* var_scale;
svd_var_scale = svd_var_scale .* var_scale;
gpca_var_scale = gpca_var_scale .* var_scale;
lpca_var_scale = lpca_var_scale .* var_scale;
epca_var_scale = transpose(epca_var_scale) .* var_scale;

%Packing scaled variance matrices into one matrix
var_scale_plot = [pca_var_scale;nnmf_var_scale;gpca_var_scale;lpca_var_scale;epca_var_scale];
%var_scale_plot = [pca_var_scale;nnmf_var_scale;svd_var_scale;gpca_var_scale;lpca_var_scale;epca_var_scale];


%Normalizing pca time and scaling other matrices by this factor
time_scale = 1 ./ pca_time;
pca_time_scale = pca_time .* time_scale;
nnmf_time_scale = nnmf_time .* time_scale;
svd_time_scale = svd_time .* time_scale;
gpca_time_scale = gpca_time .* time_scale;
lpca_time_scale = lpca_time .* time_scale;
epca_time_scale = epca_time .* time_scale;

%Packing scaled time values into one matrix
time_scale_plot = [pca_time_scale,nnmf_time_scale,gpca_time_scale,lpca_time_scale,epca_time_scale];
%time_scale_plot = [pca_time_scale,nnmf_time_scale,svd_time_scale,gpca_time_scale,lpca_time_scale,epca_time_scale];

%% Packaging into mat file

cd(dataDir);

%Creation and packing of struct with name, output, var, scale, time
pca_struct = {};
pca_struct.name_info = 'Names of algorithm in order of appearance in struct';
pca_struct.name{1} = ('PCA');
pca_struct.output_info = 'Reduced matrix (rpinciple components) from each algorithm';
pca_struct.output{1} = pca_PC;
pca_struct.var_info = 'Variance (eigenvalues) from each algorithm';
pca_struct.var{1} = pca_var;
pca_struct.name{2} = ('NMF');
pca_struct.output{2} = nnmf_PC;
pca_struct.var{2} = vecnorm(nnmf_var);
%pca_struct.name{3} = ('SVD');
%pca_struct.output{3} = svd_PC;
%pca_struct.var{3} = svd_var;
pca_struct.name{3} = ('gPCA');
pca_struct.output{3} = gpca_PC;
pca_struct.var{3} = vecnorm(gpca_var);
pca_struct.name{4} = ('lPCA');
pca_struct.output{4} = lpca_PC;
pca_struct.var{4} = vecnorm(lpca_var);
pca_struct.name{5} = ('ePCA');
pca_struct.output{5} = epca_PC;
pca_struct.var{5} = epca_var;
pca_struct.var_scale_info = 'Scaled variances (scaled by normalizing factor of PCA PC1), used for easy comparison across algos';
pca_struct.var_scale = var_scale_plot;
pca_struct.time_info = 'Run time of each algorithm';
pca_struct.time = algo_time;
pca_struct.time_scale_info = 'Scaled times (scaled by normalizing factor of PCA time), used for easy comparison across algos';
pca_struct.time_scale = time_scale_plot;

%filename = strcat('pca_process_output_', int2str(binF));

%Saving mat file
save pca_process_output pca_struct
%save(filename,'-struct','pca_struct');
end