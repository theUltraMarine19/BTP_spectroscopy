function [Y, A, S] = compute_fly_batch(X, lambda, l, psz, mu, num_patches_each_img, patch_indices, clean_indices, num_batches, D) % l for intermediate dimension

% define the number of batches and batch size
num_each_batch = floor(num_patches_each_img/num_batches);

% initialize strictly positive matrices
if (isempty(D))
    A = abs(rand(size(X,1), l));
else
    A = D;
end
S = abs(randn(l, num_patches_each_img));
nextS = abs(randn(l, num_patches_each_img));

% rescale A to have unit column norm
init_norms = sqrt(sum(A.^2,1));
A = bsxfun(@rdivide,A, init_norms);
old_A = A;
old_S = S;

if min(X(:) < 0)
    disp("Positivity of X violated")
    exit(1)
end

% Define scaling constants
scaleX = max(X(:));

% to avoid overflow/underflow, scale the matrix
X = X./max(X(:));

% construction of phi
phi_cell = cell(numel(patch_indices), 1);
% for i = 1:numel(patch_indices)
%     temp_matrix = zeros(numel(clean_indices{i}), size(X,1));
%     for j = 1:numel(clean_indices{i})
%         temp_matrix(j, clean_indices{i}(j)) = 1;
%     end
%     phi_cell{i} = temp_matrix;
% end


% % Declaring of phi*X and phi*A
% phi_A = cell(numel(patch_indices), 1);
% phi_X = cell(numel(patch_indices), 1);
% for iter = 1:size(X,2)
%     phi_A{iter} = phi_cell{iter} * A;
%     phi_X{iter} = phi_cell{iter} * X(:, iter);
% end

% % Storing of phi_A'*phi_A into a matrix
% Atranspose_A = zeros(size(A,1), size(A,1), numel(patch_indices));
% for iter = 1:size(X,2)
%     Atranspose_A(:,:,iter) = phi_A{iter}' * phi_A{iter};
% end

% for i = 1:numel(patch_indices)
%     ele_J(i) = 0.5 * sum((phi_X{i} - phi_A{i}*S(:,i)).^2);
% end
% J = sum(ele_J) + lambda*sum(S(:));
% old_J = 2*J;  %random initialization
% JB = sum(ele_J) + lambda*sum(S(:));

%% Hoyer's NNSC Algorithm
% mu_adj1 = 0;
% mu_adj2 = 0;
% mu_adj3 = 0;
red_factor = 0.25;

% multiplier_tmp = cell(numel(patch_indices),1);
tmpA = zeros(size(A));

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
    
    phi_A = cell(numel(patch_indices_tmp), 1);
    phi_X = cell(numel(patch_indices_tmp), 1);
    multiplier_tmp = cell(numel(patch_indices_tmp),1);
    ele_J = zeros(numel(patch_indices_tmp), 1);
    
    for i = 1:numel(patch_indices_tmp)
        limit = numel(clean_indices{(batch_iter - 1)*num_each_batch + i});
        temp_matrix = sparse(1:limit, clean_indices{(batch_iter - 1)*num_each_batch + i}, ones(1, limit), limit, size(X,1));
        phi_cell{(batch_iter - 1)*num_each_batch + i} = temp_matrix;
    end

    
    % initialize J's here
    parfor i = 1:numel(patch_indices_tmp)
        phi_X_tmp = phi_cell{(batch_iter - 1)*num_each_batch + i} * X_tmp(:,i);
        phi_A_tmp = phi_cell{(batch_iter - 1)*num_each_batch + i} * A;
        ele_J(i) = 0.5 * sum((phi_X_tmp - phi_A_tmp*S_tmp(:,i)).^2);
    end
    J = sum(ele_J) + lambda*sum(S_tmp(:));
    old_J = 2*J;  %random initialization

    parfor iter = 1:size(X_tmp,2)
        phi_X{iter} = phi_cell{(batch_iter - 1)*num_each_batch + iter} * X_tmp(:,iter);
    end

    while abs((J - old_J)/old_J) > 1e-4
    
        for iter = 1:size(X_tmp,2)
            phi_A{iter} = phi_cell{(batch_iter - 1)*num_each_batch + iter} * A;
        end

        multiplier = zeros(size(A));
        for iter=1:numel(patch_indices_tmp)
            multiplier_tmp{iter} = phi_cell{(batch_iter - 1)*num_each_batch + iter}' * (phi_A{iter}*S_tmp(:,iter) - phi_X{iter})*S_tmp(:,iter)';
            multiplier = multiplier + multiplier_tmp{iter};
        end

        tmpA = A - mu*multiplier;
        tmpA(tmpA < 0) = 0;
        norms = sqrt(sum(tmpA.^2,1));
        tmpA = bsxfun(@rdivide,tmpA,norms);
        old_A = A;
        A = tmpA;

        tmp1 = zeros(size(A,2), size(X_tmp,2));
        tmp2 = zeros(size(A,2), size(S_tmp,2));
        parfor iter = 1:size(X_tmp,2)
            tmp1(:, iter) = phi_A{iter}' * phi_X{iter};
            tmp2_inter = phi_A{iter}' * phi_A{iter};  % can be made FASTER
            tmp2(:, iter) = tmp2_inter * S_tmp(:, iter);
        end
        tmp2 = tmp2 + lambda;
    %     tmp2 = A'*A*S + lambda;
        nextS_tmp = S_tmp.*tmp1./tmp2;
        old_S_tmp = S_tmp;
        S_tmp = nextS_tmp;

        old_J = J;
        ele_J = zeros(numel(patch_indices_tmp), 1);
        parfor i = 1:numel(patch_indices_tmp)
            ele_J(i) = 0.5 * sum((phi_X{i} - phi_A{i}*S_tmp(:,i)).^2);
        end
        J = sum(ele_J) + lambda*sum(S_tmp(:))

        if J > old_J
            mu = mu * red_factor;
            A = old_A;
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
Y = A * S;
Y = Y.*scaleX;

% Crucial debugging
numel(find(isnan(A)));
disp("fly error")
J
