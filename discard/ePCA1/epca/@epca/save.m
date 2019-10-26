% function save(a)
%
% saves 'a' if the last save is more than 'a.saveinterval' away
% (default 30min)
%
function a = save(a,void)

% do not save just store 
if a.aux.lastsave == 0
  a.aux.lastsave = clock;
end

if etime(clock,a.aux.lastsave) > a.aux.saveinterval
  
  if nargin > 1
    tmpW = zeros(numel(void),size(a.W,2));
    tmpW(~void,:) = a.W;
    W = a.W;
    a.W = tmpW;
  end
  
  a.aux.lastsave = clock;
  save(a.savefile,'a');
  log(a,sprintf('Saved the paramaters in file ''%s'' on %s\n',a.savefile,datestr(now)),1);

  if nargin > 1
    a.W = W;
  end
  
end
