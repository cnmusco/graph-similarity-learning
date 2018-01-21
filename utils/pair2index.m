function [ edgeIndex ] = pair2index( n,u,v )
    % Given a vertex set size n and vertex indices u,v \in [n] with u < v, 
    % returns the index of the row corresponding to (u,v) in the edge vertex
    % incidence matrix, in the canonical ordering used in all our code, including genB.m

    edgeIndex = (u-1)*(n-u/2)+(v-u);
end

