function [P, PP,L] = HMM_E_step(logPxs,Trans,Pinit);

[S,T] = size(logPxs);

alpha  = zeros(S,T);     % alpha messages
gamma  = zeros(S,T);     % gamma messages
cFact  = zeros(1,T);

for t = 1:T;    % FORWARD
    fact = max(logPxs(:,t));
    lPxs = logPxs(:,t)-fact;
    if t==1;  alpha(:,t) = exp(lPxs) .* Pinit;
    else;     alpha(:,t) = exp(lPxs) .* ( Trans*alpha(:,t-1) );    end                
    cFact(t) = log( sum(alpha(:,t)) + realmin ) + fact ;
    alpha(:,t) = alpha(:,t)/sum(alpha(:,t));
end

L = sum(cFact);
gamma(:,T) = alpha(:,T); % BACKWARD
for t = (T-1):-1:1;
        PP = gamma(:,t+1) ./ ( Trans*alpha(:,t)+realmin);
        PP = PP' * Trans;
        gamma(:,t) = PP'.*alpha(:,t);
end

PP = gamma(:,2:T) ./ ( Trans*alpha(:,1:(T-1)) +realmin);
PP = PP * alpha(:,1:(T-1))';
PP = PP .* Trans;

P  = gamma;
