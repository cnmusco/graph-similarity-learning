function [ L ] = generateRandomLaplacian(n)
    % Generates a random Laplacian for an n node graph for testing purposes.
    L = rand(n,n);
    L =1/2*(L+L');
    L = (L>0.5)+0.0;
    L = L - diag(diag(L));
    s = sum(L);
    L = diag(s)-L;
end

