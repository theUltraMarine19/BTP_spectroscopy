k = 25;   % num components
psz = 8;  % patch size

%     old_image = img(1:112,1:31);
%     old_image = mat2gray(old_image)*255;
%     num_img_patches = (size(old_image,1)-psz+1)*(size(old_image,2)-psz+1);
%     % [~, X, X1, ~, ~, ~, ~] = data_gmm(old_image, psz, num_img_patches, 0, 1, 'complete', 'random');
%     [~, X, ~, ~] = data_test(old_image, psz, num_img_patches, 0);
%     [label, model, llh] = mixGaussEm(X, k);
%     for i=1:size(model.mu,2)
%         disp(det(model.Sigma(:,:,i)));
%     end
%     save(strcat('gmm-train_pureParaffin_', num2str(k), '.mat'), 'model', 'X');

    corr = [0.2 0.5 0.8];

    test_img = double(img(:,:,:));
    test_img = mat2gray(test_img)*255;

for iter = 1:3

    num_img_patches = (size(test_img,1)-psz+1)*(size(test_img,2)-psz+1)*(size(test_img,3)-psz+1);
    [new_img, new_img_interp, Y, Y1, Y_interp, patch_indices, ~, clean_indices, ~] = data_gmm(test_img, psz, num_img_patches, corr(iter), 1, 'complete', 'random');
%   [new_img, X_test, patch_indices, clean_indices] = data_test(test_img, psz, num_img_patches, corr(iter));

    % In case you want to train on interpolated image
    [label, model, llh] = mixGaussEm(Y_interp, k);
    for i=1:size(model.mu,2)
        disp(det(model.Sigma(:,:,i)));
    end
    %   save(strcat('gmm-train_pureParaffin_Fly_', num2str(k), '.mat'), 'model', 'Y_interp');
    %     
    %   load a pre-trained model
    %   load(strcat('gmm-train_pureParaffin_', num2str(k), '.mat'), 'model', 'Y_interp');

    % Adding noise    
    avg_intensity = sum(test_img(:))/numel(test_img);
    noise_var = (0.01 * avg_intensity)^2;

    % init_norms = sqrt(sum(Y.^2,1));
    % Y = bsxfun(@rdivide,Y, init_norms);

%   [Z,~,~, phi_cel] = EMfly(Y, model, patch_indices, clean_indices, noise_var);
    
    % Use gmm_estimate incase you want to just predict for offline trained images
    Z = EMfly(Y, model, patch_indices, clean_indices, noise_var);
    [recon_img, recon_img1, RMSE, zero_pxls] = reconstruct_gmm(new_img, test_img, psz, Z, patch_indices, [], [], 'complete');
%   [recon_img, recon_img1, RMSE, mask] = reconstruct(new_img, test_img, psz, Z, patch_indices, [], [], 'complete');
    RMSE1 = norm(test_img(:) - new_img_interp(:),2)/norm(test_img(:),2);
    save(strcat('Si_9_model_fly_', num2str(RMSE), '_c', num2str(corr(iter)), '_psz', num2str(psz), '.mat'), 'recon_img1', 'new_img', 'test_img', 'new_img_interp', 'RMSE', 'RMSE1', 'k');
%   imwrite(mat2gray(recon_img1), strcat('s', num2str(k1), '_7_', 'c', num2str(corr(iter)), '_psz', num2str(psz), '_RMSE', num2str(RMSE), '.png'));
end