function [resistances,u,v] = getRes(weights)
    % Gets the effective resistences of a graph given its edge weights.
    %
    % usage : 
    %
    % input :
    %
    %  * weights: (n choose 2) length vector of edge weights. 
    %      The ordering of the weights must correspond to the canonical edge
    %      ordering used in utils/pair2index.m and utils/genB.m.
    %
    % output :
    %
    %  * resistances: (n choose 2) length vector of effective resistances.
    %
    %  * u,v: (n choose 2) length vectors containing the node indices for
    %  each effective resistance pair.
    
    n = ceil(sqrt(2*size(weights,1)));
    
    L = w2L(weights);
    L = pinv(L);
    
    resistances = zeros(size(weights));
    u = zeros(size(weights));
    v = zeros(size(weights));
    for i=1:n
        for j=i+1:n
            resistances(pair2index(n,i,j)) = L(i,i)+L(j,j)-2*L(i,j);
            u(pair2index(n,i,j)) = i;
            v(pair2index(n,i,j)) = j;
        end
    end
    
end