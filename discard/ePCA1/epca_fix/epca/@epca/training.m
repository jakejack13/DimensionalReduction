function [dout,alg] =  training(a,d)


switch a.optimizer
 case 'armijo'
  [dout,alg] = training_armijo(a,d);
  return;
 case 'gd2'
  log(a,'No normalizing step\n',1);
end

X = get_x(d);
void = sum(X)==0;
valid = ~void;

nDims_old = size(X,2);
X(:,void) = [];

[nPoints,nDims] = size(X);
lastobjective = -Inf;
converged = 0;

%
% reload old params or...
% 
if a.algorithm.use_prev_train  
  if sum(sum(a.W(void,:)))
    error('not zero entries where there are no words\n');
  end
  a.W = a.W(~void,:);
  if size(a.H,2) ~= size(X,1)
    error('dimension mismatch of count data and a.H');
  end
  if size(a.W,1) ~= size(X,2)
    error('dimension mismatch of count data and a.W');
  end 
else  % ... initialize new ones 
  a.W = randn(nDims,a.dim)/(nDims);
  a.H = randn(a.dim,nPoints)/(a.dim*nPoints);
  switch a.optimizer
    case 'gd'
     norms = sqrt(sum(a.W'.^2,2));
     a.W = a.W./(ones(size(a.W,1),1) * norms');
  end
end


a.aux.lastsave = clock;
starttime = cputime;
perm = [];

last_totalobjective = -Inf;
totalobjective = 0;
t = 0;

if a.batchsize > 0
  [dout,alg] = training_newton(a,d);
  return;
end

%
%
% Main algorithm
%
%
ddH = zeros(size(a.H));
ddW = zeros(size(a.W));

if a.iterations == 0
  % start with two gradient steps for H to get the scale appr right
  %[objective,dH] = epca_obj(X,a.W,a.H,a.distr,0,1);
  [objective,dH] = gradH(a,X);
  a.H = a.H + a.etaH * dH;
end
oldW = a.W;
oldH = a.H;


%[totalobjective,dW] = epca_obj(X,a.W,a.H,a.distr,1);
[totalobjective,dW] = gradW(a,X);


lastobjective = totalobjective;
conv_count = 0;

while (a.iterations < a.maxIterations) || a.iterations > 10
  
  Wold = a.W;
  ddW = a.etaW * dW + a.momW * ddW;
  a.W = a.W + ddW;
  switch a.optimizer
    case 'gd'
     norms = sqrt(sum(a.W.^2));
     a.W = a.W./(ones(size(a.W,1),1) * norms);
  end
  %[objective,dH] = epca_obj(X,a.W,a.H,a.distr,0,1);
  [objective,dH] = gradH(a,X);
  
  conv = (objective - lastobjective)/nPoints;
  if conv < 0 || sum(isnan(a.W(:))) || sum(isinf(a.W(:))) || ...
	isnan(objective) || isnan(conv)
    a.W = Wold;
    if isnan(objective)
      %[objective,dH] = epca_obj(X,a.W,a.H,a.distr,0,1);
      [objective,dH] = gradH(a,X);
    end
    ddW = zeros(size(a.W));
    a.etaW = .8 * a.etaW;
    log(a,sprintf('change = %04g\tlowered etaW to %g\n',conv,a.etaW),2);
    conv = Inf;
    reject_w = 1;
  else 
    a.etaW = 1.01 * a.etaW;
    reject_w = 0;
    paramchange = sum(sum(abs(a.W-Wold)));
    lastobjective = objective;
    log(a,['Iter ', num2str(a.iterations), ' updW, objective function : ' num2str(lastobjective/nPoints) ', conv: ' num2str(conv) ' change in W: ' num2str(paramchange) '\n'],2);
  end

  Hold  = a.H;
  ddH = a.etaH * dH + a.momH * ddH;
  a.H = a.H + ddH;
  %[objective,dW] = epca_obj(X,a.W,a.H,a.distr,1);
  [objective,dW] = gradW(a,X);
  
  conv = (objective - lastobjective)/nPoints;
  
  if conv < 0 || sum(isnan(a.H(:))) || sum(isinf(a.H(:))) || ...
	isnan(objective) || isnan(conv)
    a.H = Hold;
    if isnan(objective)
      %[objective,dW] = epca_obj(X,a.W,a.H,a.distr,1,0);
      [objective,dW] = gradW(a,X);
    end
    ddH = zeros(size(a.H));
    a.etaH = .8 * a.etaH;
    log(a,sprintf('change = %04g\tlowered etaH to %g\n',conv,a.etaH),2);
    conv = Inf;
    reject_h = 1;
    if isnan(objective), objective=-Inf; end
  else 
    reject_h = 0;
    a.etaH = 1.01 * a.etaH;
    paramchange = sum(sum(abs(a.H-Hold)));
    lastobjective = objective;
    log(a,['Iter ', num2str(a.iterations), ' updH, objective function : ' num2str(lastobjective/nPoints) ', conv: ' num2str(conv) ' change in H= ' num2str(paramchange) '\n'],2);
  end
  

  a.iterations = a.iterations+1;

  if ~reject_w & ~reject_h
    conv = (lastobjective-totalobjective)/nPoints;
    if (conv < a.eps) && a.iterations > 10 && conv_count == 1
      log(a,sprintf('final change in objective %g\n',conv),1);
      converged = 1;
      break;
    elseif (conv < a.eps) && a.iterations > 10
      conv_count = 1;
    else
      conv_count = 0;
    end
    totalobjective = lastobjective;
  else
    conv_count = 0;
  end
  
  % ... save data so now and then (30min)
  a = save(a,void);
  
end  

log(a,sprintf('Traing time %g hours\n',(a.algorithm.training_time+cputime-starttime)/3600));

if converged
  log(a,sprintf('Converged after %d iterations, final objective: %g\n',a.iterations,objective/nPoints));
else  
  log(a,sprintf('Maximum of iterations (%d) reached\n',a.iterations));
end

alg = a;
alg.objective_value = objective/nPoints;
alg.W = zeros(numel(void),a.dim);
alg.W(valid,:) = a.W;
alg.H = a.H;
alg.iterations = a.iterations;
alg.converged = converged;

alg.algorithm.use_prev_train = 1;
alg.algorithm.training_time =alg.algorithm.training_time + cputime-starttime;

dout = data(a.H');
if numel(d.Y)
  dout.Y = d.Y;
end


function [obj,dW] = gradW(a,X);

WH=transpose(a.H)*transpose(a.W);
expWH = exp(WH);
dW=a.H*X-a.H*expWH;
obj = epca_objective(X,WH,expWH);
%obj = sum(sum(mex_sparsedotmult(X,WH)-expWH));
dW = transpose(dW);


%tic; WH=(W*H0)'; expWH=exp(WH);tt=X-expWH; ddW =H0*tt; obj = sum(sum(mex_sparsedotmult(X,WH)-expWH)); ddW=ddW';toc;
%tic; [totalobjective,dW] = epca_obj(X,a.W,a.H,'poisson',1,0);toc

function [obj,dH] = gradH(a,X);


WH=transpose(a.H)*transpose(a.W);
expWH = exp(WH);
dH=X*a.W-expWH*a.W;
%obj = sum(sum(mex_sparsedotmult(X,WH)-expWH));
obj = epca_objective(X,WH,expWH);
dH = transpose(dH);

%tic; WH=(W*H0)'; expWH=exp(WH);tt=X-expWH; ddH =tt*a.W; obj = sum(sum(mex_sparsedotmult(X,WH)-expWH)); ddW=ddW';toc;
%tic; [totalobjective,dW] = epca_obj(X,a.W,a.H,'poisson',0,1);toc
