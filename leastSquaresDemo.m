%Generate a small, clustered k-nearest neighbor graph.
n = 100;
X = .3*randn(n,3);
X(1:n/2,:) = X(1:n/2,:)+repmat([.6,0,0],n/2,1);
X(n/2+1:end,:) = X(n/2+1:end,:)+repmat([0,.6,0],n/2,1);
A = knn_sym(X,7);
L = diag(sum(A))-A;

%the true weights and effective resistances
wstar = L2w(L);
[rstar,u,v] = getRes(wstar);

% sample 2000 effective resistances to use as constraint set.
rcont = sparse(length(rstar),1);
cind = randi(length(rstar),2000,1);
rcont(cind) = rstar(cind);

% set up 4 paralell threads for running the parfor loops
myCluster = parcluster('local');
myCluster.NumWorkers = 4;
parpool(4)
% test the method for small graphs, using full set of constraints to check
% convergence
[wtildes,rtildes,optErrs,Aerrs] = effResGDSmall(rstar, L, 0, .001*ones(size(wstar)), 2000);

%test the method for larger graphs, with stochastic coordinate descent
%instead of full gradient computations. Again use full set of constraints.
[wtilde rtilde] = effResGD(rstar,u,v, L, 0, .001*ones(size(wstar)), 1000,0,3000, 'GDLS');