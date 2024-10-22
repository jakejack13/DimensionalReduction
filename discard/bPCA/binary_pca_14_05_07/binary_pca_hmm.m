function [U,V,Delta,LogLik,Q,hmm] = binary_pca_hmm(X,L,C,iid,epsil)
%
% [U,V,Delta,LogLik,Q,hmm] = binary_pca_hmm(X,L,C,iid,epsil)
% 
% Estimates a mixture of Binary PCA data model, optionally embedded in a Hidden Markov Model  
%
% Inputs:
% X (NxD) : binary data matrix, N observations of D variables, order of observations is used when HMM is estimated
% L       : dimensionality of latent space, may be zero.
% C       : number of clusters > 0
% iid     : learn mixture if iid==1, else learn HMM
% epsil   : threshold on relative loglik. increase for stopping the iterations
%
% Outputs:
% U (NxLxC)   : L dimensional latent coordinates of the N observations on the C binary PCA models
% V (LxDxC)   : LxD dimensional projection matrices of the C binary PCA models that map latent coordinates to log-odds
% Delta (CxD) : the D-dimensional off-sets of the log-odds of the C binary PCA models
% LogLik      : sequence of data log-likelihoods throughout the EM iterations
% Q (NxC)     : posterior probabilities of the N observations to be generated by each one of the C binary PCA models  
% hmm 	      : HMM transition and initial probabilities 
%
% Jakob Verbeek, 2007, INRIA Rhone-Alpes, France.
% Bin. PCA algorithm    A. Schein, L. Saul and L. Ungar. , AI and Stat's 2003
%

[N,D] = size(X); 

fprintf('--> %d data points with %d features loaded.\n',N,D);  
if iid; fprintf('--> learning mixture'); else    fprintf('--> learning HMM');end;
fprintf(' with %d mixture components\n',C);
fprintf('--> computing %d dimensional binary PCA\n',L);

x = X;
X = 2*X-1;

% initialization
Delta = randperm(N); Delta = X(Delta(1:C),:)/100; % initialize the mean vectors from random data vectors
U     = 1e-4*randn(N,L,C);                        %                factors
V     = 1e-4*randn(L,D,C);                        %            and loading matrices
for c=1:C; Th(:,:,c) = U(:,:,c)*V(:,:,c) + ones(N,1)*Delta(c,:); end; % PCA log-odds parameters

Q     = normalise(rand(N,C),2);                   % initialize the component posteriors
hmm.Pinit  = ones(C,1)/C;                         %            the HMM initial distribution
hmm.Ptrans = ones(C)/C;                           %            the HMM transition matrix

max_iters = 100;
LogLik=[];

for iter=1:max_iters; % EM loops
    
    % E-step: assignment:  component conditional likelihood
    for c=1:C; LogPxc(:,c) = sum( x.*Th(:,:,c) + log(sigmoid(-Th(:,:,c)))  , 2 );  end;
    % E-step: assignment:  combine with state prior
    [Q,PP,LogLik(iter)] = HMM_E_step(LogPxc',hmm.Ptrans,hmm.Pinit); Q=Q';
    
    % M-step: re-estimate state prior
    if iid; hmm.Pinit = mean(Q,1)'; hmm.Ptrans = repmat(hmm.Pinit,1,C);
    else    hmm.Pinit = Q(1,:)';   hmm.Ptrans = normalise(PP,1);end;
              
    % M-step: update latent projections
    T= tanh(Th/2) ./ Th;  %compute new bounds
    for c=1:C;
        B = V(:,:,c) * (X - T(:,:,c).*(ones(N,1)*Delta(c,:)))';
	for n=1:N;	% U update
            A = (V(:,:,c).* (ones(L,1)*T(n,:,c))) * V(:,:,c)';
            U(n,:,c) = mldivide(A,B(:,n))';
        end;
        Th(:,:,c)    = U(:,:,c)*V(:,:,c) + ones(N,1)*Delta(c,:);             % PCA log-odds parameters
    end;
    
    % M-step: update loading matrices and mixing weights
    T= tanh(Th/2) ./ Th;  %compute new bounds
    for c=1:C;
        U2 = [U(:,:,c) ones(N,1)]; U2 = U2.*repmat(Q(:,c),1,L+1);
        B  = U2' * X;
        for d=1:D;	% V update 
            A = (U2'.* (ones(L+1,1)*T(:,d,c)')) * [U(:,:,c) ones(N,1)];        
            V2 = mldivide(A,B(:,d));
            Delta(c,d) = V2(end);
            if L>0; V(:,d,c)   = V2(1:end-1);end;
        end;
        Th(:,:,c)    = U(:,:,c)*V(:,:,c) + ones(N,1)*Delta(c,:);             % PCA log-odds parameters
    end;
    
    if iter>1; rel_incr = ( LogLik(end)-LogLik(end-1) ) / abs(LogLik(end-1));else rel_incr=1;end;
    fprintf('%3d rel_incr: %f  log-lik: %f\n',iter,rel_incr,LogLik(end)/numel(X));
    if iter>2 && rel_incr < epsil; break; end;
    
end; % EM iterations

% post-processing : de-correlate latent variables and set projection vectors to unit length 
for c=1:C;
	Delta(c,:) = Delta(c,:) + mean(U(:,:,c))*V(:,:,c);
	U(:,:,c)   = U(:,:,c) - ones(N,1)*mean(U(:,:,c));
	T = U(:,:,c)*V(:,:,c);
    [u,s,v] = svds(T,L);
    U(:,:,c) = u * s;	V(:,:,c) = v' ;
end;
