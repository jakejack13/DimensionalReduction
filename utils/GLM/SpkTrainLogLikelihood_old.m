function L = SpkTrainLogLikelihood_old(q,f)

% L = SpkTrainValuation_old(q,f)
% 
% computes log-likelihood of spike train S with intensity function f.
% INPUTS:
%     q: binned spike train
%     f: predicted intensity function (in number of spikes)
% 
% OUTPUT:
%     L: log-likelihood
%
% Adrien Peyrache, 2014 (following Harris, 2004)
%
% changed to log2() - Luke Sjulson, 2016

L = -f+q.*log2(f);
L = sum(L(f>0));
if isinf(L) || isnan(L)
   L =  0;
end
