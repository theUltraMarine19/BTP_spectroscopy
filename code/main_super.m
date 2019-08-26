% for super-resolution via interpolation
function X = main_super(img, scene)
lambda = 2e-4;
mu = 3e-5;
num_patches_each_img = 46 * 46;
corrupt_pattern = 3;
psz = 6;
l = 12;

%% Main Script

% generate the regularly corrupted image, and the subsmapled 'knowns' image
[new_img, X, X1, X_super, patch_indices, patch_indices_unfilled, clean_indices, clean_indices_unfilled, super_indices] = data_test_hyper(img, psz, num_patches_each_img, corrupt_pattern, 'complete', 'regular', 'interp');

% generate initial guess dictionary on X_super
[Y_super, A_init, S_init] = compute_batch_noPhi(X_super, lambda, l, mu, numel(super_indices), super_indices, 4);

% use this to learn the main dictionary
[Y, A, S] = compute_fly_batch(X, lambda, l, psz, mu, num_patches_each_img, patch_indices, clean_indices, 4, A_init);

[recon_img2, recon_img21, RMSE2, zero_pxls2] = reconstruct(new_img, img, psz, Y_super, patch_indices, [], patch_indices_unfilled, 'partial');
[recon_img, recon_img1, RMSE, zero_pxls] = reconstruct(new_img, img, psz, Y, patch_indices, [], patch_indices_unfilled, 'partial');

RMSE
save(strcat(scene, '_RMSE_', num2str(RMSE), '_c', num2str(corrupt_pattern), '_psz', num2str(psz), '_onTheFly', '.mat'), 'recon_img1', 'new_img', 'img', 'RMSE2', 'zero_pxls', 'mu', 'lambda', 'l', 'num_patches_each_img', 'A', 'A_init');
