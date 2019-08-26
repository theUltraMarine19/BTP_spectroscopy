function [Y, A, S] = compute_fly_batch_poisson(X, lambda, l, psz, mu, num_patches_each_img, patch_indices, clean_indices, num_batches) % l for intermediate dimension

% define the number of batches and batch size
num_each_batch = floor(num_patches_each_img/num_batches);

% initialize strictly positive matrices
A = abs(rand(size(X,1), l));
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

%% Hoyer's NNSC Algorithm
red_factor = 0.25;

% multiplier_tmp = cell(numel(patch_indices),1);
tmpA = zeros(size(A));

for batch_iter=1:num_batches
    
    disp(strcat('batch ', num2str(batch_iter), ' is starting')); 
    
    %% Division of batches
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
    
    %% Declare auxiliary variables
    phi_A = cell(numel(patch_indices_tmp), 1);
    phi_X = cell(numel(patch_indices_tmp), 1);
    multiplier_tmp = cell(numel(patch_indices_tmp),1);
    
    %% Construction of phi's
    for i = 1:numel(patch_indices_tmp)
        limit = numel(clean_indices{(batch_iter - 1)*num_each_batch + i});
        temp_matrix = sparse(1:limit, clean_indices{(batch_iter - 1)*num_each_batch + i}, ones(1, limit), limit, size(X,1));
        phi_cell{(batch_iter - 1)*num_each_batch + i} = temp_matrix;
    end

    
    %% initialize J's here
    ele_J = zeros(numel(patch_indices_tmp), 1);
    for i = 1:numel(patch_indices_tmp)
        phi_X_tmp = phi_cell{(batch_iter - 1)*num_each_batch + i} * X_tmp(:,i);
        phi_A_tmp = phi_cell{(batch_iter - 1)*num_each_batch + i} * A;
        phi_A_S_tmp = phi_A_tmp * S_tmp(:,i);
        ele_J_col = phi_X_tmp.*log((phi_X_tmp+1e-11)./phi_A_S_tmp) - phi_X_tmp + phi_A_S_tmp;
%         if numel(find(isnan(ele_J_col))) > 0
% %             disp(numel(find(isnan(phi_X_tmp.*log(phi_X_tmp./phi_A_S_tmp)))));
%             f = find(isnan(phi_X_tmp.*log(phi_X_tmp./phi_A_S_tmp)));
%             phi_X_tmp_vec = phi_X_tmp(:);
%             log_tmp = log(phi_X_tmp./phi_A_S_tmp);
%             log_vec = log_tmp(:);
%             phi_X_tmp_vec(f);
%             log_vec(f);
%             disp("ERROR");
%             numel(find(isnan(A)))
%         end
        ele_J(i) = sum(ele_J_col(:));
%         for j=1:size(phi_cell{(batch_iter - 1)*num_each_batch + i},1)
%             ele_J(i) = ele_J(i) + phi_X_tmp(j,i)*log(phi_X_tmp(j,i)/phi_A_S_tmp(j)) - phi_X_tmp(j,i) + phi_A_S_tmp(j);
%         end
    end
    J = sum(ele_J) + lambda*sum(S_tmp(:))
    old_J = 2*J;  %random initialization
    
    %% Construct phi*X
    for iter = 1:size(X_tmp,2)
        phi_X{iter} = phi_cell{(batch_iter - 1)*num_each_batch + iter} * X_tmp(:,iter);
    end
    
    %% Gradient descent
    while abs((J - old_J)/old_J) > 2.5e-4
        
%         A(1:10, 1:10)
        
        % Update phi*A
        for iter = 1:size(X_tmp,2)
            phi_A{iter} = phi_cell{(batch_iter - 1)*num_each_batch + iter} * A;
        end
        
        %% update of A
        multiplier = zeros(size(A));
        for iter=1:numel(patch_indices_tmp)
            multiplier_tmp{iter} = phi_cell{(batch_iter - 1)*num_each_batch + iter}' * (1 - phi_X{iter}./(phi_A{iter}*S_tmp(:,iter)))*S_tmp(:,iter)';
%             cnt = numel(find(phi_A{iter}*S_tmp(:,iter) == 0))
%             if cnt ~= 0
%                 cnt
%             end
            multiplier = multiplier + multiplier_tmp{iter};
        end
        
%         multiplier(1:10, 1:10)
        tmpA = A - mu*multiplier;
        
        % Projection part of gradient descent
        tmpA(tmpA < 0) = 0;
        norms = sqrt(sum(tmpA.^2,1));
        tmpA = bsxfun(@rdivide,tmpA,norms);
        old_A = A;
        A = tmpA;
        
        %% Update of S
        tmp1 = zeros(size(A,2), size(X_tmp,2));
        tmp2 = zeros(size(A,2), size(S_tmp,2));
        for iter = 1:size(X_tmp,2)
            tmp1(:, iter) = phi_A{iter}' * (phi_X{iter}./(phi_A{iter} * S_tmp(:, iter)));
            tmp2(:, iter) = phi_A{iter}' * ones(size(phi_cell{(batch_iter - 1)*num_each_batch + iter},1),1);
        end
        tmp2 = tmp2 + lambda;
        nextS_tmp = S_tmp.*tmp1./tmp2;
        old_S_tmp = S_tmp;
        S_tmp = nextS_tmp;
        
        %% Update of cost
        old_J = J;
        ele_J = zeros(numel(patch_indices_tmp), 1);
        for i = 1:numel(patch_indices_tmp)
            phi_A_S_tmp = phi_A{i} * S_tmp(:,i);
            ele_J_col = phi_X{i}.*log((phi_X{i}+1e-11)./phi_A_S_tmp) - phi_X{i} + phi_A_S_tmp;
%             if numel(find(isnan(ele_J_col))) > 0
% %                 disp(numel(find(isnan(phi_X{i}.*log(phi_X{i}./phi_A_S_tmp)))));
%                 f = find(isnan(phi_X{i}.*log(phi_X{i}./phi_A_S_tmp)));
%                 phi_X_tmp_vec = phi_X{i}(:);
%                 log_tmp = log(phi_X{i}./phi_A_S_tmp);
%                 log_vec = log_tmp(:);
%                 phi_X_tmp_vec(f(1));
%                 log_vec(f(1));
%                 phi_A_S_tmp(f(1));
%                 disp("ERROR");
%                 numel(find(isnan(A)))
%             end
            ele_J(i) = sum(ele_J_col(:));
        end
        J = sum(ele_J) + lambda*sum(S_tmp(:))
        
        %% Check for overshooting of local minima
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
