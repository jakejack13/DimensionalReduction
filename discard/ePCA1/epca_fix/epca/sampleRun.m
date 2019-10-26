
% For 1-dim this algorithm should find the largest eigenvector of
% the data
X = [.5* randn(100,1),1.5*randn(100,1)];
a = epca({'distr=''normal''','dim=1'});
a.minLogL = 1e-4;
[d,e] = train(a,data(sparse(X)));

e.A % solution of the epca algorithm
[U,D] = eig(X'*X);
U(:,2) % largest eigenvector


clear
load 20News_w100
d = data(Xtest',ytest');
a = epca({'distr=''poisson''','dim=10','optimizer=''smd'''});
% lower the learning rates, the parameters are not so stable as
% exp(A*V) is taken during computation
a.etaW = 0.05;
a.etaH = 0.05;
a.minLogL = 1e-4;
[dout,e] = train(a,d);

clf;
plot(e,[0:.1:6])
hold on;
scatter3(X(:,1),X(:,2),X(:,3))
