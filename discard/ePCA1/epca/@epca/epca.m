function a = epca(hyper)

%=========================================================================
% Exponential Family PCA spider object
%=========================================================================
% a = epca(hyperParam)
% 
% Generates an EPCA object with given hyperparamters.
%
% Hyperparameters (with defaults)
%  k=1                  -- number of neighbours
%  batch=1              -- true if to be computed in batch 
%                           (requires more memory)
%  child=kernel;        -- kernel function (distances computation/ 
%                           default is linear kernel)
%  output_preimage=0    -- whether output index from training sample 
%                           of preimage 
%                          instead of actual label
%
% Model
%  dat                  -- data used for neighbours computation
%
% Methods:
%  train, test
%
% XXX example
%========================================================================
% Reference : A Generalization of Principal Component Analysis to
%             the Exponential Family
% Author    : M.Collins, S.Dasgupta, R.E.Schapire
% Link      : 
%========================================================================

% hyperparams

a.dim = 10;
a.H = [];
a.W = [];
a.maxIterations = 1e+4;
a.distr = 'normal'; %'bernoulli','poisson'

a.iterations = 0;

a.etaW = .002;
a.etaH = .002;
a.momW = .8;
a.momH = .8;

a.mu = 0.001;
a.lambdaW = 0.99;
a.lambdaH = 0.7;

a.eps = 0.0001; 

a.normalize = 0;
a.retrain = 0;

a.optimizer = 'smd'; 

a.batchsize = 0;
a.converged = 0;
a.objective_value = -Inf;

a.savefile = strrep([num2str(clock) num2str(rand)],' ','');
a.aux.lastsave = 0; % will cause saveing 1/2 hour after learning started
a.aux.saveinterval = 1800; % save every 1/2 hour

a.aux.dW = [];
a.aux.vW = [];
p=algorithm('epca');
a= class(a,'epca',p);

if nargin==1,
    eval_hyper;
end;

