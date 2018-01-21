% generates a full (n choose 2) x n vertex edge incidence matrix for a graph
% with n vertices. Edges are in the canonical order defined by pair2index.m.
function [ B ] = genB(n)
    m = n*(n-1)/2;
    B = sparse(m,n);
    for i=1:n
         for j=i+1:n
             B(pair2index(n,i,j),i) = 1;
             B(pair2index(n,i,j),j) = -1;
         end
    end
end

