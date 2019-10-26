% [ind,perm,t] = nextbatch(a,ind,N)
% 
% retrieves the next a.batchsize indices from ind and reshuffles
% the indices of the dataset if necessary. In case a.batchsize == 0
% the whole dataset will be the batch. 
%
% t   1 if a pass through the data is complete, 0 otherwise
%
function [ind,perm,t] = nextbatch(a,perm,N)

t = 0;
bs = a.batchsize;

if bs==0 || bs >= N;
  ind = 1:N;
  perm = [];
  t = 1;
  return;
end

if numel(perm) == 0
  perm(end+1:end+N) = randperm(N);
% reshuffle dataset and append
elseif numel(perm) <= bs
  t = 1; 
  perm(end+1:end+N) = randperm(N);
end

% retrieve the next batch
ind = perm(1:bs);

% and remove it from the rest
perm(1:bs) = [];


