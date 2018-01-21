function [ M ] = w2A( w )
    % Convert an (n choose 2) length vector of edge weights to an associated
    % symmetric adjacency matrix.
    %
    % Note that the ordering of the input vector matters. It is consistent with
    % the ordering used pair2index.m and genB.m.

    n = ceil(sqrt(2*size(w,1)));
    M = zeros(n,n);
    for i=1:n
        for j=i+1:n
            M(i,j) = w(pair2index(n,i,j));
            M(j,i) = w(pair2index(n,i,j));
        end
    end
end
