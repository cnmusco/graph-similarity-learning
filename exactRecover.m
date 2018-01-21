function [ w ] = exactRecover( r,lambda )
    % Given a full set of (n choose 2) effective resitances, r, recovers the edge
    % weights of the unique underlying graph with resistances matching r. 
    %
    % Implements the formula in Theorem 1 of "Learning Networks from Random
    % Walk-Based Node Similarities".
    %
    % If r is a valid set of effective resistances, set the regularization
    % parameter lambda = 0 to recover the exact underlying graph. If r is
    % estimated or noisy, use positive lambda to heuritically generate a set 
    % of (possibly negative) weights which may be thresholded to generate a
    % graph with effective resistances approximately matching r. 

    % usage : 
    %
    % input :
    %
    %  * r: (n choose 2) length vector of pairwise effective resistances. 
    %       Note that the ordering of the input vector matters. It should be
    %       consistent with the ordering used by pair2index.m and e.g. A2w.m
    %
    %  * lambda: regularization parameter, lambda >= 0.
    %
    % output :
    %
    %  * w: (n choose 2) length vector of edge weights. If r is a valid set of
    %  effective resistances, w(i) >= 0 for all i, up to roundoff error.

    R = w2A(r);
    n = size(R,1);
    ONES = 1/sqrt(n)*ones(n,1);
    Linv = R - ONES*(ONES'*R);
    Linv = Linv - (Linv * ONES)*ONES';
    L = inv(-Linv/2 + ONES * ONES' + lambda*eye(n,n));
    L = L - (L*ONES)*ONES';
    L = L - ONES*(ONES'*L);
    w = L2w(L);
end