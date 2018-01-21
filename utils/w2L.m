function [ L ] = w2L( w )
    % Convert an (n choose 2) length vector of edge weights to a corresponding 
    % graph Laplacian
    %
    % Note that the ordering of the input vector matters. It is consistent with
    % the ordering used pair2index.m and genB.m.

    A = w2A(w);
    L = diag(sum(A))-A;
end

