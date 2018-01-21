function [ w ] = L2w( L )
%Convert a Laplacian matrix to its associated edge weights.
    n = size(L,1);
    m = n*(n-1)/2;
    w = zeros(m,1);
    for i=1:n
        for j=i+1:n
            w(pair2index(n,i,j)) = -L(i,j);
        end
    end
end

