function a = knn_sym(x,k)
% generates a symmetric k nearest neighbors graph from a set of points x.
%

    [n,~] = size(x);

    k = min(k,n-1);

    list = zeros(k,n);

    for i = 1:n
      y = x - ones(n,1)*x(i,:);
      s = sum(y.^2,2);  
      [val, ord] = sort(s);
      if (val(2) > eps) % no repeat point
        list(:,i) = ord(2:(k+1));
      else
        [~,start] = min(val < eps);
        list(:,i) = ord([start:(start+k-1)]);
      end
    end

    a = sparse(list(:), kron([1:n]', ones(k,1)), 1, n, n);

    a = (a | a')*1;

end

