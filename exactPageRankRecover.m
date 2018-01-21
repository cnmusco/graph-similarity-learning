function [ w ] = exactPageRankRecover( P,alpha,lambda )
    % Given an n x n matrix P of pairwise personalized PageRank scores, computed 
    % with jump back probability alpha, recovers the edge weights of the unique
    % underlying graph with PageRank scores matching P. 
    %
    % Implements the formula in Theorem 3 of "Learning Networks from Random
    % Walk-Based Node Similarities".
    %
    % If P is a valid set of pairwise PageRank scores, set the regularization
    % parameter lambda = 0 to recover the exact underlying graph. If P is
    % estimated or noisy, use positive lambda to heuritically generate a set 
    % of (possibly negative) weights which may be thresholded to generate a
    % graph with pairwise PageRank scores approximately matching P. 

    % usage : 
    %
    % input :
    %
    %  * P: n x n matrix of pairwise personalized PageRank scores.
    %
    %  * alpha: jump back probability that the scores in P were computed with.
    %
    %  * lambda: regularization parameter, lambda >= 0.
    %
    % output :
    %
    %  * w: (n choose 2) length vector of edge weights. If P is a valid set of
    %  personalized PageRank scores, w(i) >= 0 for all i, up to roundoff error.

    n = size(P,1);
    AD = 1/(1-alpha)*((1+alpha)*eye(n,n) - 2*alpha*inv(P+lambda*eye(n,n)));
    
    for j=1:n
        AD(:,j) = AD(:,j)/max(AD(:,j));
    end
    w = L2w(-AD);
end


