% to fill in missing pixels (assurance that all pixels are reconstructed)
function [Y,S] = test(A, X, lambda, patch_indices, clean_indices) % is phi to be externally supplied (derive from patches of the image)
S = abs(randn(size(A,2), size(X,2)));
old_S = S;

psz = sqrt(size(X, 1));
scaleX = max(X(:));

% to avoid overflow/underflow, scale the matrix
X = X./scaleX;

% construction of phi
phi_cell = cell(numel(patch_indices), 1);
for i = 1:numel(patch_indices)
    temp_matrix = zeros(numel(clean_indices{i}), psz*psz);
    for j = 1:numel(clean_indices{i})
        temp_matrix(j, clean_indices{i}(j)) = 1;
    end
    phi_cell{i} = temp_matrix;
end

ele_J = zeros(numel(patch_indices), 1);
for i = 1:numel(patch_indices)
    ele_J(i) = 0.5 * sum((phi_cell{i}*X(:,i) - phi_cell{i}*A*S(:,i)).^2);
end
J = sum(ele_J) + lambda*sum(S(:));
old_J = 2*J;  %random initializatiom

while abs((J - old_J)/old_J) > 1e-6
%     tmp1 = A'*X;    % TROUBLE in this step
    tmp1 = zeros(size(A,2), size(X,2));
    tmp2 = zeros(size(A,2), size(S,2));
    for iter = 1:size(X,2)
        tmp_mult1 = phi_cell{iter}*A;
        tmp_mult1_transpose = tmp_mult1';
        tmp_mult2 = phi_cell{iter}*X(:, iter);
        tmp1(:, iter) = tmp_mult1_transpose * tmp_mult2;
        tmp2_inter = tmp_mult1_transpose * tmp_mult1;
        tmp2(:, iter) = tmp2_inter * S(:, iter);
    end
%     tmp2 = A'*A*S + lambda;
    nextS = S.*tmp1./tmp2;
    S = nextS;
    old_J = J;
    ele_J = zeros(numel(patch_indices), 1);
    for i = 1:numel(patch_indices)
        ele_J(i) = 0.5 * sum((phi_cell{i}*X(:,i) - phi_cell{i}*A*S(:,i)).^2);
    end
    J = sum(ele_J) + lambda*sum(S(:))
end

disp("test");
J

% obtain Y = AS
Y = A * S;
Y = Y.*scaleX;

end
        
        

