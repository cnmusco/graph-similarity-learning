function [P] = getPageRank(weights,alpha)
    % Computes all personalized Pagerank scores of a graph given its edge weights.
    %
    % usage : 
    %
    % input :
    %
    %  * weights: (n choose 2) length vector of edge weights. 
    %      The ordering of the weights must correspond to the canonical edge
    %      ordering used in utils/pair2index.m and utils/genB.m.
    %
    %  * alpha: the jump pack probability used in personalized PageRank computation
    %
    % output :
    %
    %  * P: n x n matrix contained all pairwise personalized PageRank scores.
    %
    n = ceil(sqrt(2*size(weights,1)));
    
    A = w2A(weights);
    Dinv = diag(1./sum(A));
    % the random walk matrix
    W = .5*(eye(n,n) + A*Dinv);
    P = alpha*inv(eye(n,n)-(1-alpha)*W);    
end