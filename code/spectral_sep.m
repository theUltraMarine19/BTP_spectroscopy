function [Y, D1, S] = spectral_sep(X, lambda, l, psz, mu, num_patches_each_img, patch_indices, clean_indices, num_batches, D2) % l for intermediate dimension
% D2 is the known dictionary that is supplied

% define the number of batches and batch size
num_each_batch = floor(num_patches_each_img/num_batches);

% initialize strictly positive matrices
D1 = abs(rand(size(X,1), l));
S = abs(randn(2*l, num_patches_each_img));
nextS = abs(randn(2*l, num_patches_each_img));

% rescale A to have unit column norm
init_norms = sqrt(sum(D1.^2,1));
D1 = bsxfun(@rdivide,D1, init_norms);
old_D1 = D1;
old_S = S;

if min(X(:) < 0)
    disp("Positivity of X violated")
    return;
end

% Define scaling constants
scaleX = max(X(:));

% to avoid overflow/underflow, scale the matrix
X = X./max(X(:));


%% Hoyer's NNSC Algorithm
% mu_adj1 = 0;
% mu_adj2 = 0;
% mu_adj3 = 0;
red_factor = 0.25;

% multiplier_tmp = cell(numel(patch_indices),1);
tmpD1 = zeros(size(D1));

for batch_iter=1:num_batches
    
    disp(strcat('batch ', num2str(batch_iter), ' is starting')); 
    % Division of batches
    if batch_iter~=num_batches
        X_tmp = X(:,(batch_iter - 1)*num_each_batch+1:batch_iter*num_each_batch);
        S_tmp = S(:,(batch_iter - 1)*num_each_batch+1:batch_iter*num_each_batch);
        nextS_tmp = nextS(:,(batch_iter - 1)*num_each_batch+1:batch_iter*num_each_batch);
        old_S_tmp = old_S(:,(batch_iter - 1)*num_each_batch+1:batch_iter*num_each_batch);
        patch_indices_tmp = patch_indices(:,(batch_iter - 1)*num_each_batch+1:batch_iter*num_each_batch);
    else
        X_tmp = X(:,(batch_iter - 1)*num_each_batch+1:size(X,2));
        S_tmp = S(:,(batch_iter - 1)*num_each_batch+1:size(S,2));
        nextS_tmp = nextS(:,(batch_iter - 1)*num_each_batch+1:size(nextS,2));
        old_S_tmp = old_S(:,(batch_iter - 1)*num_each_batch+1:size(old_S,2));
        patch_indices_tmp = patch_indices(:,(batch_iter - 1)*num_each_batch+1:size(old_S,2));
    end
    
    
    % initialize J's here
    J = 0.5 * sum(sum((X_tmp - D1 * S_tmp(1:l,:) - D2 * S_tmp((l+1):2*l,:)).^2)) + lambda*sum(S_tmp(:));
    old_J = 2*J;  %random initialization

    
    while abs((J - old_J)/old_J) > 1e-4
    
        multiplier = (D1*S_tmp(1:l,:) + D2*S_tmp((l+1):2*l,:) - X_tmp)*S_tmp(1:l,:)';
        tmpD1 = D1 - mu*multiplier;
        tmpD1(tmpD1 < 0) = 0;
        norms = sqrt(sum(tmpD1.^2,1));
        tmpD1 = bsxfun(@rdivide,tmpD1,norms);
        old_D1 = D1;
        D1 = tmpD1;

        A = horzcat(D1, D2);
        tmp1 = A'*X_tmp;
        tmp2 = A'*A*S_tmp + lambda;
        nextS_tmp = S_tmp.*tmp1./tmp2;
        old_S_tmp = S_tmp;
        S_tmp = nextS_tmp;

        old_J = J;
        J = 0.5 * sum(sum((X_tmp - D1 * S_tmp(1:l,:) - D2 * S_tmp((l+1):2*l,:)).^2)) + lambda*sum(S_tmp(:))

        if J > old_J
            mu = mu * red_factor;
            D1 = old_D1;
            S_tmp = old_S_tmp;
        end

    end
    if batch_iter~=num_batches
       S(:,(batch_iter - 1)*num_each_batch+1:batch_iter*num_each_batch) = S_tmp;
    else
       S(:,(batch_iter - 1)*num_each_batch+1:size(S,2)) = S_tmp;
       
    end
    
end
% obtain Y = AS
Y = D1 * S(1:l, :);
Y = Y.*scaleX;

% Crucial debugging
numel(find(isnan(A)));
disp("fly error")
J
