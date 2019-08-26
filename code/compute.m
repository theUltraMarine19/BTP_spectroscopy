function [A, S] = compute(X, lambda, l, mu, psz, num_patches_each_img) % l for intermediate dimension (To be tuned)

% initialize strictly positive matrices
A = abs(rand(size(X,1), l)); %*(max(X(:)) - min(X(:))) + min(X(:));
% S = abs(randn(l, 20*10*num_patches_each_img));
S = abs(randn(l, num_patches_each_img));

% rescale A to have unit column norm
init_norms = sqrt(sum(A.^2,1));
A = bsxfun(@rdivide,A, init_norms);
old_A = A;
old_S = S;

if min(X(:) < 0)
    disp("Positivity of X violated")
    exit(1)
end

% to avoid overflow/underflow, scale the matrix
X = X./max(X(:));

J = 0.5 * sum(sum((X - A*S).^2)) + lambda*sum(S(:));
old_J = 2*J;  %random initializatiom

%% code for individual update of A, S
% update S for 100 steps
% for t = 1:100
%     tmp1 = A'*X;
%     tmp2 = A'*A*S + lambda;
%     nextS = S.*tmp1./tmp2;
%     S = nextS;
% end

% update A for 100 steps
% for t = 1:100
%     multiplier = (A*S - X)*S';
%     tmpA = A - mu*multiplier;
%     tmpA(tmpA < 0) = 0;
%     norms = sqrt(sum(tmpA.^2,1));
%     tmpA = bsxfun(@rdivide,tmpA,norms);
%     A = tmpA;
% end

%% Hoyer's NNSC Algorithm
% mu_adj1 = 0;
% mu_adj2 = 0;
% mu_adj3 = 0;
red_factor = 0.5;

while abs((J - old_J)/old_J) > 5e-4
%     disp("Running NNSC...");
    
    multiplier = (A*S - X)*S';
    tmpA = A - mu*multiplier;
    tmpA(tmpA < 0) = 0;
    norms = sqrt(sum(tmpA.^2,1));
    tmpA = bsxfun(@rdivide,tmpA,norms);
    old_A = A;
    A = tmpA;
    
    tmp1 = A'*X;
    tmp2 = A'*A*S + lambda;
    nextS = S.*tmp1./tmp2;
    old_S = S;
    S = nextS;
    
    old_J = J;
    J = 0.5 * sum(sum((X - A*S).^2)) + lambda*sum(S(:))
    
    % Variable step-size
%     if abs((J - old_J)/old_J) < 1e-3 && mu_adj1 == 0
%         mu = mu/2;
%         mu_adj1 = 1;
%     elseif abs((J - old_J)/old_J) < 1e-4 && mu_adj2 == 0
%         mu_adj2 = 1;
%         mu = mu/5;
%     elseif abs((J - old_J)/old_J) < 1e-5 && mu_adj3 == 0
%         mu = mu/2;
%         mu_adj3 = 1;
%     end

    if J > old_J
        mu = mu * red_factor;
        A = old_A;
        S = old_S;
    end
    
    
end

% Crucial debugging
numel(find(isnan(A)))
disp("train")
J

% Visualization 
% for iter=1:size(A,2)
%     figure;
%     imshow(mat2gray(reshape(A(:,iter), [psz psz])));
% end


