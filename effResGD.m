function [wtildes,rtildes,optErrs,Aerrs] = effResGD(r, u,v, Ltrue, lambda, w0, iters, degreeConst, batchSize, method )
% Given a set of effective resistance constraints and regularization
% parameter lambda, uses gradient descent to try to find a graph which
% minimizes the squared norm difference from these constraints plus lambda*tr(L).

% Input:
%   r: effective resistance target vector. (n choose 2) length sparse
%   vector with all nonzeros considered as constraints.

%   u,v: (n choose 2) length vectors containing ids of edge end points in 
%        canonical edge ordering. These vectors can be obtained by the
%        utils/getRes.m method and are always the same for a given n. Passed in for efficiency.

%   Ltrue: For evaluation, you may pass in the true Laplacian to measure
%   error to it. If the true Laplacian is unknown, just pass in zeros(n,n).

%   lambda: regularization parameter (will minimize ||rtilde -r||_2^2 + lambda*tr(L)
%           default = 0.

%   w0: starting weights (default all ones)

%   iters:  list of iteration counts, [iter1, iter2, ...]. will output the
%           solution at each of iter1, iter2, ... (default [100,500,1000])
%           for evaluation purposes

%   degreeConst: if this is 1 it will in each iteration restrict each node to
%                have the same degree that it  has in Ltrue. default is 0.
%                NOTE: CURRENTLY UNIMPLEMENTED. Setting degreeConst = 1
%                will have no effect.

%   batchSize: For large graphs, it may be neccesary to use coordinate
%              descent. batchSize = number of edges sampled and updated at 
%              each step. Default = m, which causes full gradient computations

%   method: 'GD' for gradient descent, 'GDLS' for gradient descent with
%            line search (default)

% Output:
%   wtildes: The learned edge weights, saved at each iteration count in iters.

%   rtildes: The effective resistances of the learned graph, saved at each
%            iteration count in iters.

%   optErrs: The objection function value, recorded at each iteration count
%            in iters

%   AErrs:   The distance between the learned adjacency matrix and the true
%            adjacency matrix (if known) in Frobenius norm, saved at each
%            iteration count in iters.

m = length(u);
n = ceil(sqrt(2*m));

if nargin < 10
    method = 'GDLS';
end
if nargin < 9
    batchSize = m;
end
if nargin  < 8
    degreeConst  = 0;
end
if nargin < 7
    iters = [100,500,1000];
end
% initial complete graph
if nargin < 6
    w0 = ones(m);
end
if nargin < 5
    lambda = 0;
end

%The constrained edges. Any edge with r(i) = 0, we ignore.
cind = find(r~=0);
numConstraints = length(cind);

%will store the current effective resistances at the constrained edges
rtilde = zeros(numConstraints,1);
%will store the current errors from the given effective resistances
deltaw = zeros(numConstraints,1);

wtilde = w0;
Ltilde = w2L(w0);

%%
%output variables.
wtildes = zeros(length(iters),length(wtilde));
rtildes = zeros(length(iters),length(rtilde));
%%
% instrumentation variables
% measures the objective function value -- the distance from the true resistances
optErrs = zeros(1,max(iters));
% measures the distance between the recovered adjacency  matrix and the true
% adjacency matrix (if input for testing purposes) in Frobrenius norm. 
Aerrs = zeros(1,max(iters));
Atrue = diag(diag(Ltrue)) - Ltrue;

linesearch = strcmp(method,'GDLS');
%shift for the line search starting point.
shift  = 0;
numSearchs = 10;
if(strcmp(method,'GD') || linesearch)
        %the indexs of the edges to update in each iteration. By default we
        %update all unless the batchSize is set. In which case we update a
        %random subset.
        R = zeros(batchSize,numConstraints);
        %special case when batchSize = m we will perform deterministic full
        %gradient computations.
        wInd = (1:m)';
        for i = 1:max(iters)
            % choose random batch of edges to update
            if(batchSize < m)
                wInd = randsample(1:m,batchSize);
            end
 
            Linv = inv(Ltilde+1);
            parfor j=1:numConstraints
                edgeInd = cind(j);
                rtilde(j) = Linv(u(edgeInd),u(edgeInd))+Linv(v(edgeInd),v(edgeInd))-2*Linv(u(edgeInd),v(edgeInd));
            end
            
            % the current errors for the constained edges.
            deltaw = r(cind) - rtilde;  
            if(mod(i,10)==0)
                display(['On iteration: ',num2str(i), ' error is ', num2str(norm(deltaw)/norm(r(cind)))]);
                fprintf('Number of entries in wtilde above 0: %d\n',sum(wtilde > 0));
            end
            %run gradient step
            for k = 1:numConstraints
                colVec = Linv(:,u(cind(k)))-Linv(:,v(cind(k)));
                for j = 1:length(wInd)
                    R(j,k) = colVec(u(wInd(j)))-colVec(v(wInd(j)));
                end
            end
            grad = sparse(wInd,[1],2*R.^2*deltaw + 2*lambda,m,1); 
            % line search for step size. 
            % TODO: maybe unroll loop for efficiency
            if(linesearch)
                optVals = zeros(numSearchs,1);
                for j = 1:numSearchs
                    wopt = wtilde - grad/(norm(grad)*1.5^(2*(j+shift)-numSearchs));
                    wopt = (wopt > 0).* wopt;
                    Ltilde = w2L(wopt);
                    Linv = inv(Ltilde+1);
                    rtildeOpt = zeros(numConstraints,1);
                    for k=1:numConstraints
                        edgeInd = cind(k);
                        rtildeOpt(k) = Linv(u(edgeInd),u(edgeInd))+Linv(v(edgeInd),v(edgeInd))-2*Linv(u(edgeInd),v(edgeInd));
                    end
                    optVals(j) = norm(r(cind)-rtildeOpt)^2  +2*lambda*sum(wopt);
                end
                [optmin, optj] = min(optVals);
                %Uncomment for diagnostics:
                %fprintf('The best iteration of line search was: %d\n',optj');
                optErrs(i) = optmin./norm(r(cind))^2;
                wnew = wtilde - grad/(norm(grad)*1.5^(2*(optj+shift)-numSearchs));
                if(optj == numSearchs)
                    shift = shift+5;
                    fprintf('The line search shift is now: %d\n',shift');
                elseif(optj == 1)
                    shift = shift-5;
                    fprintf('The line search shift is now: %d\n',shift');
                end 
            else
                wnew = wtilde - 20*grad/iter;
            end
            wtilde = (wnew>0).*wnew;
            Atilde = w2A(wtilde);
            Ltilde  = diag(sum(Atilde))-Atilde;
            wtilde = L2w(Ltilde);
            
            if(ismember(i,iters))
                wtildes(find(iters == i),:) = wtilde;
                rtildes(find(iters == i),:) = rtilde;
                save(sprintf('iter%d',i),'wtildes','rtildes','Aerrs','optErrs');
            end
            Aerrs(i) = norm(Atilde - Atrue, 'fro')/norm(Atrue, 'fro');
            if(mod(i,10)==0)
                display(['On iteration: ',num2str(i), ' A err is ', num2str(Aerrs(i))]);
            end
        end
else
    error('GDHeuristic:BadInput','the specificed method was not recognized')
end
Linv = inv(w2L(wtilde)+1);
for k=1:numConstraints
    edgeInd = cind(k);
    rtilde(k) = Linv(u(edgeInd),u(edgeInd))+Linv(v(edgeInd),v(edgeInd))-2*Linv(u(edgeInd),v(edgeInd));
end

end

