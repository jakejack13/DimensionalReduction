function makeSpeedModel_GLM_v21(SSDbasedir)

% function makeSpeedModel_GLM_v21(SSDbasedir)
%
% another script to fit GLM models. This uses a 25 ms binsize and incorporates
% a time lag between HPC and NAc. Also uses all significant SVDs from SVD_v1.mat.
%
% Luke Sjulson, 2016-12-31




% % for testing
% clear all
% close all
% tic
% SSDbasedir = '/home/luke/Dropbox/forAdrien/data/BF1-1_CPP_2014-10-08';

%% start of actual function
%Parameters:
binT        =  0.025; % length of time bins for binning spikes
NAcLag      =  0.030; % best fit as per Adrien's testing
rateThresh  =  0.01; % cells below this firing rate will be discarded (entries will be NaNs)
Nvalid        = 10; % number of cross-validations to perform (split data into this many chunks)

% load everything in
[~, basename] = fileparts(SSDbasedir);
cd(SSDbasedir);
load BehavEpochs.mat;
load GeneralInfo_v1.mat;
load pos.mat;
load NAcSpikes.mat;
load ThreeChamberStats_v2.mat;
load SVD_v1.mat

if binT~=SVD_v1.binT, error('bin sizes do not match'); end


%% bin and smooth

% relevant intervalsets
cocEp = ThreeChamberStats_v2.cocaineInt;
salEp = ThreeChamberStats_v2.salineInt;
bothEp     = union(cocEp, salEp);

% Restrict everything to the two chambers
% NAcSpikes.S = Restrict(NAcSpikes.S, bothEp);
linSpd      = Restrict(pos.linSpd, bothEp);
% posTs       = median(diff(Range(pos.linSpd)));


%% Spike time binning
Q       = MakeQfromS(NAcSpikes.S, binT);
Q = Restrict(Q, bothEp);
dQ      = Data(Q);


%% make location indicator vectors

% 1 in saline location, 0 otherwise
% dSloc = zeros(length(Range(Q)), 1);
% [~, ix] = Restrict(Q, salEp);
% dSloc(ix) = 1;

% 1 in cocaine location, 0 otherwise
dCloc = zeros(length(Range(Q)), 1);
[~, ix] = Restrict(Q, cocEp);
dCloc(ix) = 1;

%% Interpolating speed
dSpd     = nanzscore(interp1(Range(linSpd), Data(linSpd), Range(Q), 'linear', 'extrap'));

%% interpolating PCs (with timelag)
nPCs = SVD_v1.HPYRmanualEigs;
dPCs = zeros(size(dQ, 1), nPCs);

for idx = 1:nPCs
   dPCs(:,idx) = interp1(SVD_v1.timevec+NAcLag, SVD_v1.uHPYR(:,idx), Range(Q), 'nearest');
end

if isnan(dPCs)
   warning(sprintf('Interpolated PCs for %s had NaNs...will set to zero.', basename));
   dPCs(isnan(dPCs)) = 0;
end

%% doing cross-validation

beta58 = NaN(1+nPCs, length(NAcSpikes.S), Nvalid);
beta59 = NaN(2+nPCs, length(NAcSpikes.S), Nvalid);
beta60 = NaN(2+nPCs, length(NAcSpikes.S), Nvalid);
beta61 = NaN(3+nPCs, length(NAcSpikes.S), Nvalid);

beta58Pval = NaN(1+nPCs, length(NAcSpikes.S), Nvalid);
beta59Pval = NaN(2+nPCs, length(NAcSpikes.S), Nvalid);
beta60Pval = NaN(2+nPCs, length(NAcSpikes.S), Nvalid);
beta61Pval = NaN(3+nPCs, length(NAcSpikes.S), Nvalid);

L58 = NaN(Nvalid, length(NAcSpikes.S));
L59 = NaN(Nvalid, length(NAcSpikes.S));
L60 = NaN(Nvalid, length(NAcSpikes.S));
L61 = NaN(Nvalid, length(NAcSpikes.S));


% split intervals up into Nvalid chunks, loop over each
setNum = ceil(rand(length(Range(Q)), 1)*Nvalid);

for crossIdx = 1:Nvalid
   
   % extract the data for the training intervals
   dQTrain = dQ(setNum~=crossIdx, :);
   %   dSlocTrain = dSloc(setNum~=crossIdx, :);
   dClocTrain = dCloc(setNum~=crossIdx, :);
   dSpdTrain = dSpd(setNum~=crossIdx, :);
   dPCTrain = dPCs(setNum~=crossIdx, :);
   
   
   % data for the test intervals
   dQTest = dQ(setNum==crossIdx, :);
   %   dSlocTest = dSloc(setNum==crossIdx, :);
   dClocTest = dCloc(setNum==crossIdx, :);
   dSpdTest = dSpd(setNum==crossIdx, :);
   dPCTest = dPCs(setNum==crossIdx, :);
   
   
   % loop over all cells
   for cellIdx = 1:length(NAcSpikes.S)
      % for cellIdx = 1:5 % for testing only
      
      % only fit model if cell's firing rate is high enough
      if mean(dQ(:,cellIdx)./binT) > rateThresh
         
         
         % model 58: rate + PCs
         % log10(y) = B(1) + B(2)*PC1 + ... + B(n+1)*PCn
         [tempBeta, ~, stats] = glmfit(dPCTrain, dQTrain(:, cellIdx), 'poisson');
         beta58(:, cellIdx, crossIdx) = tempBeta;
         beta58Pval(:, cellIdx, crossIdx) = stats.p;
         N = size(dPCTest,1);
         f58 = exp(tempBeta(1)) .* exp(sum(repmat(tempBeta(2:end)', N, 1).*dPCTest, 2));
         L58(crossIdx, cellIdx) = SpkTrainLogLikelihood_old(dQTest(:, cellIdx), f58);
                
         % model 59: rate + speed + PCs
         % log10(y) = B(1) + B(2)*dSpd + B(3)*PC1 + ... + B(n+2)*PCn
         [tempBeta, ~, stats] = glmfit([dSpdTrain dPCTrain], dQTrain(:, cellIdx), 'poisson');
         beta59(:, cellIdx, crossIdx) = tempBeta;
         beta59Pval(:, cellIdx, crossIdx) = stats.p;
         f59 = exp(tempBeta(1) + tempBeta(2).*dSpdTest) .* exp(sum(repmat(tempBeta(3:end)', N, 1).*dPCTest, 2));
         L59(crossIdx, cellIdx) = SpkTrainLogLikelihood_old(dQTest(:, cellIdx), f59);
         
         % model 60: rate + location + PCs
         % log10(y) = B(1) + B(2)*dCloc + B(3)*PC1 + ... + B(n+2)*PCn
         [tempBeta, ~, stats] = glmfit([dClocTrain dPCTrain], dQTrain(:, cellIdx), 'poisson');
         beta60(:, cellIdx, crossIdx) = tempBeta;
         beta60Pval(:, cellIdx, crossIdx) = stats.p;
         f60 = exp(tempBeta(1) + tempBeta(2).*dClocTest) .* exp(sum(repmat(tempBeta(3:end)', N, 1).*dPCTest, 2));
         L60(crossIdx, cellIdx) = SpkTrainLogLikelihood_old(dQTest(:, cellIdx), f60);
         
         % model 61: rate + location + speed + PCs
         % log10(y) = B(1) + B(2)*dCloc + B(3)*dSpd + B(4)*PC1 + ... + B(n+3)*PCn
         [tempBeta, ~, stats] = glmfit([dClocTrain dSpdTrain dPCTrain], dQTrain(:, cellIdx), 'poisson');
         beta61(:, cellIdx, crossIdx) = tempBeta;
         beta61Pval(:, cellIdx, crossIdx) = stats.p;
         f61 = exp(tempBeta(1) + tempBeta(2).*dClocTest + tempBeta(3).*dSpdTest) .* exp(sum(repmat(tempBeta(4:end)', N, 1).*dPCTest, 2));
         L61(crossIdx, cellIdx) = SpkTrainLogLikelihood_old(dQTest(:, cellIdx), f61);
         
      end
   end
end


%% copy to struct
SG.SSDbasedir = SSDbasedir;
SG.info = 'Testing multiple GLM models against each other using log likelihood, smaller bin size, NAc lag';
SG.binT = binT; % time bins in seconds
SG.NAcLag = NAcLag; % NAc spikes are delayed by this much
SG.NAcLag_info = 'NAc spikes shifted by this much to account for HPC->NAc delay';
SG.rateThresh = rateThresh;  % cells below this firing rate will be discarded (entries will be NaNs)
SG.beta_info = 'betaN is for model N, matrix is nBetas x nCells x nCrossValidations';

SG.beta58 = beta58;
SG.beta58Pval = beta58Pval;
SG.L58 = L58;
SG.beta59 = beta59;
SG.beta59Pval = beta59Pval;
SG.L59 = L59;
SG.beta60 = beta60;
SG.beta60Pval = beta60Pval;
SG.L60 = L60;
SG.beta61 = beta61;
SG.beta61Pval = beta61Pval;
SG.L61 = L61;

SpeedModel_GLM_v21 = SG;

cd(SSDbasedir);

save('SpeedModel_GLM_v21.mat', 'SpeedModel_GLM_v21');



% toc
