
clear
load 20news_w500_counts;
ind = randperm(size(counts,2));
nTrain = .5;
nDocs = size(counts,2);
Xtrain = counts(:,ind(1:ceil(nTrain*nDocs)));
ytrain = labels(ind(1:ceil(nTrain*nDocs)));
Xtest = counts(:,ind(ceil(nTrain*nDocs)+1:end));
ytrain = labels(ind(ceil(nTrain*nDocs)+1:end));

d = data(Xtrain',ytrain');
a = epca({'distr=''poisson''','dim=50'});
% lower the learning rates, the parameters are not so stable as
% exp(A*V) is taken during computation
a.etaW = 0.0005;
a.etaH = 0.0005;
a.eps = .005;
a.verbosity = 2;
a.optimizer = 'gd2';
[dout,e] = train(a,d);

dtest = test(e,data(Xtest',ytest'));

