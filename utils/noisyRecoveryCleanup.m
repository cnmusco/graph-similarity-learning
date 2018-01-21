function [ w ] = noisyRecoveryCleanup( wtilde, p )
    %cleans up a noisyly recovered weight vector with possibly negative
    %entries just by taking  the top p fraction of positive edges. E.g., set p = .2.

    w = wtilde .* (wtilde > 0);
    % take top p fraction of edges
    t = sort(w); t = t(floor((1-p)*length(w)));
    w = w.*(w > t);
    w = w ./ mean(w(find(w)));

end