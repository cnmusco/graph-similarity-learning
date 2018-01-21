function [ w ] = recoverMissing( r, lambda )
    % Heuristic method for recovering a graph from an incomplete set of
    % effective resistances. See Section 4.2  of "Learning Networks from Random
    % Walk-Based Node Similarities".
    %
    % usage : 
    %
    % input :
    %
    %  * r: (n choose 2) length vector of pairwise effective resistances.
    %        r(i) = 0 if the i^th effective resistance is missing.
    %
    %  * lambda: regularizatoin parameter with which to run exact recover
    %    procedure. lambda >= 0.
    %
    % output :
    %
    %  * w: (n choose 2) length vector of edge weights. Some may be
    %  negative.  Consider cleaning up with utils/noisyRecoveryCleanup.m.
    
    R = sparse(w2A(r));
    % fill in in missing effective resistances with shortest path distances.
    D = graphallshortestpaths(R);
    R = R + (R == 0).*D;
    % then just attempt exact recovery.
    w = exactRecover(A2w(R),lambda);
end