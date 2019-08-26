% to fill in missing pixels (assurance that all pixels are reconstructed)
function [Y,S] = test_fast_poisson(A, X, lambda, patch_indices, clean_indices) % is phi to be externally supplied (derive from patches of the image)
S = abs(randn(size(A,2), size(X,2)));
old_S = S;

scaleX = max(X(:));

% to avoid overflow/underflow, scale the matrix
X = X./scaleX;

% construction of phi
phi_cell = cell(numel(patch_indices), 1);
for i = 1:numel(patch_indices)
    limit = numel(clean_indices{i});
    temp_matrix = sparse(1:limit, clean_indices{i}, ones(1, limit), limit, size(X,1));
    phi_cell{i} = temp_matrix;
end

% Caching of phi*X and phi*A
phi_A = cell(numel(patch_indices), 1);
phi_X = cell(numel(patch_indices), 1);
for iter = 1:size(X,2)
    phi_A{iter} = phi_cell{iter} * A;
    phi_X{iter} = phi_cell{iter} * X(:, iter);
end

ele_J = zeros(numel(patch_indices), 1);
for i = 1:numel(patch_indices)
    phi_A_S_tmp = phi_A{i} * S(:,i);
    ele_J_col = phi_X{i}.*log((phi_X{i}+1e-11)./phi_A_S_tmp) - phi_X{i} + phi_A_S_tmp;
    ele_J(i) = sum(ele_J_col(:));
end
J = sum(ele_J) + lambda*sum(S(:));
old_J = 2*J;  %random initializatiom

while abs((J - old_J)/old_J) > 2.5e-4
%     tmp1 = A'*X;    
    tmp1 = zeros(size(A,2), size(X,2));
    tmp2 = zeros(size(A,2), size(S,2));
    for iter = 1:size(X,2)
        tmp1(:, iter) = phi_A{iter}' * (phi_X{iter}./(phi_A{iter} * S(:, iter)));
        tmp2(:, iter) = phi_A{iter}' * ones(size(phi_cell{iter},1),1);    
    end
    tmp2 = tmp2 + lambda;
%     tmp2 = A'*A*S + lambda;
    nextS = S.*tmp1./tmp2;
    old_S = S;
    S = nextS;
    
    old_J = J;
    ele_J = zeros(numel(patch_indices), 1);
    for i = 1:numel(patch_indices)
        phi_A_S_tmp = phi_A{i} * S(:,i);
        ele_J_col = phi_X{i}.*log((phi_X{i}+1e-11)./phi_A_S_tmp) - phi_X{i} + phi_A_S_tmp;
        ele_J(i) = sum(ele_J_col(:));
    end
    J = sum(ele_J) + lambda*sum(S(:))
    
    if J > old_J
        disp("Increase in J");
    end

end

disp("test");
J

% obtain Y = AS
Y = A * S;
Y = Y.*scaleX;

end