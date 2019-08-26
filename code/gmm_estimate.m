function Z = gmm_estimate(Y, model, signal_dim, patch_indices, clean_indices, noise_var)
    num_gmm = size(model.mu, 2)
    Z = zeros(signal_dim, size(Y,2));
    
    % construction of phi
    phi_cell = cell(numel(patch_indices), 1);
    for i = 1:numel(patch_indices)
        temp_matrix = zeros(numel(clean_indices{i}), signal_dim);
        for j = 1:numel(clean_indices{i})
            temp_matrix(j, clean_indices{i}(j)) = 1;
        end
        phi_cell{i} = temp_matrix;
    end
    
    Sigma = model.Sigma;
    mu = model.mu;
    
    for i1=1:size(Y,2)
        estimate = zeros(signal_dim,num_gmm);
        j_vals = zeros(num_gmm,1);
        phi = phi_cell{i1};
        patch_vec = phi*Y(:,i1);
        parfor i=1:num_gmm
            if det(Sigma(:,:,i)) == 0 ||  isinf(det(Sigma(:,:,i)))
                j_vals(i) = inf;
                disp("Worng place");
                continue;
            end
            estimate(:,i) = Sigma(:,:,i)*phi'*((phi*Sigma(:,:,i)*phi' + noise_var*eye(size(phi,1)))\(patch_vec - phi*mu(:,i))) + mu(:,i);
            j_vals(i) = (patch_vec - phi*estimate(:,i))'*((noise_var*eye(size(phi,1)))\(patch_vec - phi*estimate(:,i))) + (estimate(:,i) - mu(:,i))'*(Sigma(:,:,i)\(estimate(:,i) - mu(:,i))) + log(det(Sigma(:,:,i)));
        end
        [~,j_est] = min(j_vals);
        Z(:,i1) = estimate(:,j_est);
        disp([j_est i1]);
%         j_vals'
    end
end