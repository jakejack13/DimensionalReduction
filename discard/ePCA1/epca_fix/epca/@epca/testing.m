function d =  testing(a,d)

%
% fold in H keep W fixed
%
X = get_x(d);
[nPoints,nDims] = size(X);

%d.X = X*a.W;
%return;

if numel(get_y(d))
  Y = get_y(d);
  if size(Y,1) == size(X,1) & size(Y,2)==a.dim
    H = Y';
    fprintf('reusing initialization\n');
  end
end

if ~exist('H','var')
  %H = randn(a.dim,nPoints)/(a.dim*nPoints);
  H = randn(a.dim,nPoints);%/(a.dim*nPoints);
  H = H./max(abs(H(:)));
  H = max(abs(a.H(:))) * H;
end

conv = Inf;
conv_count = 0;

%[lastobjective,dH] = epca_obj(X,a.W,H,a.distr,0,1);
[lastobjective,dH] = gradH(a,H,X);

a.aux.lastsave = clock;

iter = 0;
ddH = zeros(size(dH));

while iter < a.maxIterations || iter < 10
  iter = iter + 1;
  
  Hold = H;
  dHold = dH;
  if iter < 10 % first 10 iterations to get the scale of H appr right.
    etaH =  max(H(:))/max(dH(:));
  end

  ddH = etaH * dH + a.momH * ddH;
  H = H + ddH;

  %[objective,dH] = epca_obj(X,a.W,H,a.distr,0,1);
  [objective,dH] = gradH(a,H,X);

  % the following adaption of the stepsize does not give
  % convergence guarantee but works faster in practice than the
  % Armijo rule (see testing_armijo.m)
  while (objective < lastobjective) || isnan(objective)
    log(a,sprintf('change = %g, lowered etaH to %g\n',(objective-lastobjective)/nPoints,etaH),1);
    H = Hold;
    ddH = zeros(size(dH));
    etaH = .8 * etaH;
    H = H + etaH * dHold;
    %[objective,dH] = epca_obj(X,a.W,H,a.distr,0,1);
    [objective,dH] = gradH(a,H,X);
  end
  etaH = 1.01 * etaH;

  conv = (objective - lastobjective)/nPoints;
  lastobjective = objective;
  paramchange = sum(sum(abs(Hold-H)));
  log(a,sprintf('Iter %d, objective function: %g, conv: %g, change in H : %g\n',iter,objective/nPoints,conv,paramchange));

  if conv < a.eps && conv_count == 1
    break;
  elseif conv < a.eps
    conv_count = 1;
  else
    conv_count = 0;
  end
  
  if mod(iter,100) == 0
      log(a,sprintf('Folding in on %s\n',datestr(now)),1);
  end
 
  if numel(strfind(a.savefile,'foldin')) >0
    a = save(a);
  end
  
end

log(a,sprintf('done folding in, final objective %g\n',objective/nPoints));

d.X = H';
d=set_name(d,[get_name(d) ' -> ' get_name(a)]);



function [obj,dH] = gradH(a,H,X);

WH=transpose(H)*transpose(a.W);
expWH = exp(WH);
dH=X*a.W-expWH*a.W;
%obj = sum(sum(mex_sparsedotmult(X,WH)-expWH));
obj = epca_objective(X,WH,expWH);
dH = transpose(dH);

%[objective,dH] = epca_obj(X,a.W,H,a.distr,0,1);
