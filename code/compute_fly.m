function [Y, A, S] = compute_fly(X, lambda, l, psz, mu, num_patches_each_img, patch_indices, clean_indices) % l for intermediate dimension

% initialize strictly positive matrices
A = abs(rand(size(X,1), l));
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

% Define scaling constants
scaleX = max(X(:));

% to avoid overflow/underflow, scale the matrix
X = X./max(X(:));

% construction of phi
phi_cell = cell(numel(patch_indices), 1);
for i = 1:numel(patch_indices)
    temp_matrix = zeros(numel(clean_indices{i}), size(X,1));
    for j = 1:numel(clean_indices{i})
        temp_matrix(j, clean_indices{i}(j)) = 1;
    end
    phi_cell{i} = temp_matrix;
end

% Declaring of phi*X and phi*A
phi_A = cell(numel(patch_indices), 1);
phi_X = cell(numel(patch_indices), 1);
for iter = 1:size(X,2)
    phi_A{iter} = phi_cell{iter} * A;
    phi_X{iter} = phi_cell{iter} * X(:, iter);
end

% % Storing of phi_A'*phi_A into a matrix
% Atranspose_A = zeros(size(A,1), size(A,1), numel(patch_indices));
% for iter = 1:size(X,2)
%     Atranspose_A(:,:,iter) = phi_A{iter}' * phi_A{iter};
% end

ele_J = zeros(numel(patch_indices), 1);
for i = 1:numel(patch_indices)
    ele_J(i) = 0.5 * sum((phi_X{i} - phi_A{i}*S(:,i)).^2);
end
J = sum(ele_J) + lambda*sum(S(:));
old_J = 2*J;  %random initialization

%% Hoyer's NNSC Algorithm
% mu_adj1 = 0;
% mu_adj2 = 0;
% mu_adj3 = 0;
red_factor = 0.25;

multiplier_tmp = cell(numel(patch_indices),1);
tmpA = zeros(size(A));

while abs((J - old_J)/old_J) > 3e-4
%     disp("Running NNSC...");
    for iter = 1:size(X,2)
        phi_A{iter} = phi_cell{iter} * A;
    end

%     multiplier = (A*S - X)*S';
    multiplier = zeros(size(A));
    for iter=1:numel(patch_indices)
        multiplier_tmp{iter} = phi_cell{iter}' * (phi_A{iter}*S(:,iter) - phi_X{iter})*S(:,iter)';
        multiplier = multiplier + multiplier_tmp{iter};
    end

    tmpA = A - mu*multiplier;
     tmpA(tmpA < 0) = 0;
%     tmpA(1:5, 1:5)
    norms = sqrt(sum(tmpA.^2,1));
    tmpA = bsxfun(@rdivide,tmpA,norms);
    old_A = A;
    A = tmpA;

 %     tmp1 = A'*X;    
    tmp1 = zeros(size(A,2), size(X,2));
    tmp2 = zeros(size(A,2), size(S,2));
    for iter = 1:size(X,2)
        tmp1(:, iter) = phi_A{iter}' * phi_X{iter};
        tmp2_inter = phi_A{iter}' * phi_A{iter};  % can be made FASTER
        tmp2(:, iter) = tmp2_inter * S(:, iter);
    end
    tmp2 = tmp2 + lambda;
%     tmp2 = A'*A*S + lambda;
    nextS = S.*tmp1./tmp2;
    old_S = S;
    S = nextS;

    old_J = J;
    ele_J = zeros(numel(patch_indices), 1);
    for i = 1:numel(patch_indices)
        ele_J(i) = 0.5 * sum((phi_X{i} - phi_A{i}*S(:,i)).^2);
    end
    J = sum(ele_J) + lambda*sum(S(:))

    if J > old_J
        mu = mu * red_factor;
        A = old_A;
        S = old_S;
    end

end
    
% obtain Y = AS
Y = A * S;
Y = Y.*scaleX;

% Crucial debugging
numel(find(isnan(A)));
disp("fly error")
J
