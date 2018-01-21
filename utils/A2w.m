function [ w ] = A2w( A )
%Convert an adjacency matrix to its associated edge weights.
    w = L2w(-A);
end