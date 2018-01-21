function L = sdpalgo(r)
    % Given a set of (n choose 2) effective resitances constraints, r, 
    % finds the minimum total degree graph with all effective resistances
    % upper bounded by these constraints. See Section 4.3 of ""Learning
    % Networks from Random Walk-Based Node Similarities".
    %
    % usage : 
    %
    % input :
    %
    %  * r: (n choose 2) length vector of pairwise effective resistances.
    %       If r(i) == 0, the i^th effective resistance is considered
    %       unconstrained.
    %
    % output :
    %
    %  * L: n x n graph Laplacian with minimum total degree satisfying the 
    %       effective resistance constraints

m = length(r);
n = ceil(sqrt(2*m));
B = genB(n);

c = 1;
for i = 1:length(r)
    if(r(i) ~= 0)
        measurements(c,:) = [find(B(i,:) ~= 0), r(i)];
        c = c+1;
    end
end

[len ign] = size(measurements);

cvx_begin sdp    
    variable L(n,n) symmetric semidefinite
    minimize( trace(L) )
    subject to:
	for i = 1  : len
        u = measurements(i,1);
        v =  measurements(i,2);
        ruv = measurements(i,3);
        chi = zeros(n,1);
        chi(u) = 1; 
        chi(v)=-1; 
        [L chi; chi' ruv] >=0; 
    end
cvx_end