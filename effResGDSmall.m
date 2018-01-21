function [wtildes,rtildes,optErrs,Aerrs] = effResGDSmall(r, Ltrue, lambda, w0, iters, degreeConst, method )
% Given a set of effective resistance constraints and regularization
% parameter lambda, uses gradient descent to try to find a graph which
% minimizes the squared norm difference from these constraints plus lambda*tr(L).

% Input:
%   r: effective resistance target vector. (n choose 2) length sparse
%   vector with all nonzeros considered as constraints.

%   Ltrue: For evaluation, you may pass in the true Laplacian to measure
%   error to it. If the true Laplacian is unknown, just pass in zeros(n,n).

%   lambda: regularization parameter (will minimize ||rtilde -r||_2^2 + lambda*tr(L)
%           default = 0.

%   w0: starting weights (default all ones)

%   iters:  list of iteration counts, [iter1, iter2, ...]. will output the
%           solution at each of iter1, iter2, ... (default [100,500,1000])
%           for evaluation purposes

%   degreeConst: if this is 1 it will in each iteration restrict each node to
%                have the same degree that it  has in Ltrue. default is 0.
%                NOTE: CURRENTLY UNIMPLEMENTED. Setting degreeConst = 1
%                will have no effect.

%   method: 'GD' for gradient descent, 'GDLS' for gradient descent with
%            line search (default)

% Output:
%   wtildes: The learned edge weights, saved at each iteration count in iters.

%   rtildes: The effective resistances of the learned graph, saved at each
%            iteration count in iters.

%   optErrs: The objection function value, recorded at each iteration count
%            in iters

%   AErrs:   The distance between the learned adjacency matrix and the true
%            adjacency matrix (if known) in Frobenius norm, saved at each
%            iteration count in iters.

if nargin < 7
    method = 'GDLS';
end
if nargin  < 6
    degreeConst  = 0;
end
if nargin < 5
    iters = [100,500,1000];
end
% initial complete graph
if nargin < 4
    w0 = ones(size(r));
end
if nargin < 3
    lambda = 0;
end

m = length(r);
n = ceil(sqrt(2*m));

%The constrained edges. Any edge with r(i) = 0, we ignore.
cind = find(r~=0);
constraints = length(cind);

%will store the current effective resistances at the constrained edges
rtilde = zeros(constraints,1);
%will store the current errors from the true effective resistances
deltaw = zeros(constraints,1);
B = genB(n);

wtilde = w0;
Ltilde = w2L(w0);

wtildes = zeros(length(iters),length(wtilde));
rtildes = zeros(length(iters),length(rtilde));
optErrs = zeros(1,max(iters));
Aerrs = zeros(1,max(iters));
Atrue = diag(diag(Ltrue)) - Ltrue;

numSearchs = 10;
linesearch = strcmp(method,'GDLS');
%shift for the line search starting point.
shift  = 0;
if(strcmp(method,'GD') || linesearch)
        for i = 1:max(iters)
            Linv = inv(Ltilde+1);
            rtilde = sum(B(cind,:)'.*(Linv*B(cind,:)'))';
            
            % the current errors for the constained edges.
            deltaw = r(cind) - rtilde;
            if(mod(i,50)==0)
                display(['On iteration: ',num2str(i), ' error is ', num2str(norm(deltaw)/norm(r(cind)))]);
            end
            %run gradient step
            R = B*(Linv*B(cind,:)');
            grad = 2*R.^2*deltaw + 2*lambda; 
            if(linesearch)
                optVals = zeros(numSearchs,1);
                parfor j = 1:numSearchs
                    wopt = wtilde - grad/(norm(grad)*1.5^(2*(j+shift)-numSearchs));
                    wopt = (wopt > 0).* wopt;
                    Ltilde = w2L(wopt);
                    rtildeOpt = sum(B(cind,:)'.*(inv(Ltilde+1)*B(cind,:)'))';
                    optVals(j) = norm(r(cind)-rtildeOpt)^2  +2*lambda*sum(wopt);
                end
                [optmin, optj] = min(optVals);
                optErrs(i) = optmin./norm(r(cind))^2 ;
                %Uncomment for diagnostics:
                %fprintf('The best iteration of line search was: %d\n',optj');
                wnew = wtilde - grad/(norm(grad)*1.5^(2*(optj+shift)-numSearchs));
                if(optj == numSearchs)
                    shift = shift+5;
                    fprintf('The line search shift is now: %d\n',shift');
                elseif(optj == 1)
                    shift = shift-5;
                    fprintf('The line search shift is now: %d\n',shift');
                end 
            else
                wnew = wtilde - 20*grad/iter;
            end
            wtilde = (wnew>0).*wnew;
            Atilde = w2A(wtilde);
            Ltilde  = diag(sum(Atilde))-Atilde;
            wtilde = L2w(Ltilde);
            
            if(ismember(i,iters))
                wtildes(find(iters == i),:) = wtilde;
                rtildes(find(iters == i),:) = rtilde;
                save(sprintf('iter%d',i),'wtildes','rtildes','Aerrs','optErrs');
            end
            Aerrs(i) = norm(Atilde - Atrue, 'fro')/norm(Atrue, 'fro');
            if(mod(i,50)==0)
                display(['On iteration: ',num2str(i), ' A err is ', num2str(Aerrs(i))]);
            end
        end
else
    error('GDHeuristic:BadInput','the specificed method was not recognized')
end

end

