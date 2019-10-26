
function plot(a,interv)

if numel(a.V) == 0
  error('object seems to be untrained');
end

[dim,inputDim] = size(a.W);

if inputDim > 3 
  error('can only plot up to three dimensions');
end


if nargin < 2
  interv = [0:.1:1]';
elseif size(interv,2) > size(interv,1)
    interv = interv';
end

for i=1:a.dim

  switch a.distr
   case 'normal'
    gxV = interv * a.W(i,:);
   case 'poisson'
    gxV = exp(interv * a.W(i,:));
   case 'bernoulli'
    expxv = exp(interv * a.W(i,:));
    gxV = expxv./(1+expxv);
  end
  
  switch(inputDim)
   case 1
    plot(x,gxV);
   case 2
    plot(gxV(:,1),gxV(:,2));
   case 3
    plot3(gxV(:,1),gxV(:,2),gxV(:,3),'-');
  end
  hold on
end

hold off
grid on