%% Parameters for functions

%Going from gPCA subdirectory to data directory to access data
[baseDir,dataDir] = gPCA_checkuser;
cd(dataDir) 
load('NAcSpikes.mat')

%Setting up input settings
%U = [3,2,1; 2,-3 0; 0,0,0];
%S = [5,0,0;4,0,0;2,0,0];
%t = 0:0.01:4*pi;
%V = [sin(t) + 1; sin(1.7*t) + 1; sin(5.3*t) + 1];
%I = U*S*V;

%input = round(rand(100,6));

input = extract_data(NAcSpikes,20); %binning data from raw

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

%Algorithm execution
[A,B,C] = pca(input); %PCA
Z = B * transpose(A);
[W,F] = nnmf(input, k); %NMF
Y = W * F;
[~,pU,PCs] = gPCA(input, ncol, timevec, k, M, max_iters, normalize, family); %gPCA
P = PCs * transpose(pU);
[G,E,~,D,~,~,~,~,~,~] = exp_fam_pca(input,family,k,0,k); %ePCA
Dd = D*transpose(G);
    
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

[~,lU,LCs] = lPCA(input, ncol, timevec, k, M, max_iters, normalize, family); %lPCA
L = LCs * transpose(lU);

%Moving back to base directory 
cd(baseDir);
%% Scaling tests
C_tot = 0;
for idx = 1:numel(C)
    C_tot = C_tot + C(idx);
end
C_scale = C/C_tot;

W_var = vecnorm(W);
W_tot = 0;
for idx = 1:numel(W_var)
    W_tot = W_tot + W_var(idx);
end
W_scale = W_var/W_tot;

P_var = vecnorm(PCs);
P_tot = 0;
for idx = 1:numel(W_var)
    P_tot = P_tot + P_var(idx);
end
P_scale = P_var/P_tot;

L_var = vecnorm(LCs);
L_tot = 0;
for idx = 1:numel(L_var)
    L_tot = L_tot + L_var(idx);
end
L_scale = L_var/L_tot;

E_tot = 0;
for idx = 1:numel(E)
    E_tot = E_tot + E(idx);
end
E_scale = E/E_tot;

scale = 1 / C_scale(1);
C_scale = transpose(C_scale) * scale;
W_scale = W_scale * scale;
P_scale = P_scale * scale;
L_scale = L_scale * scale;
E_scale = transpose(E_scale) * scale;

cat = categorical({'PCA','NMF','gPCA','lPCA','ePCA'});
scale_plot1 = [C_scale;W_scale;P_scale;L_scale;E_scale];
bar(cat,scale_plot1,'stacked');
%% Plotting results
close all;
f = figure;

%timeLim = 500:600;
%plotLim = 2:3;

subplot(2,2,1);
%plot(B(timeLim,plotLim));
plot(A);
title('PCA');
subplot(2,2,2);
%plot(W(timeLim,plotLim));
plot(F);
title('NMF');
subplot(2,2,3);
%plot(PCs(timeLim,plotLim));
plot(pU);
title('gPCA');
subplot(2,2,4);
%plot(D(timeLim,plotLim));
plot(G);
title('ePCA');
%plot(LCs(timeLim,plotLim));
%plot(lU);
%title('lPCA');