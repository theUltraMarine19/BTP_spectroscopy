function [] = main_fly(img, scene)
lambda = 2e-4;
mu = 3e-5;
num_patches_each_img = [5000, 8000, 10500];
corrupt_percent = [0.2 0.5 0.8];
psz = 8;
l = 12;

%% Adding noise here
% img1 = poissrnd(img, 'poisson');

%% Main Script

for iter = 1:1
    [new_img, X, X1, X_super, patch_indices, patch_indices_unfilled, clean_indices, clean_indices_unfilled, super_indices] = data_test_hyper(img, psz, num_patches_each_img(iter), corrupt_percent(iter), 'complete', 'random', 'DC');
    [Y, A, S] = compute_fly_batch(X, lambda, l, psz, mu, num_patches_each_img(iter), patch_indices, clean_indices, 2, []);
    [Y1, S1] = test_fast(A, X1, lambda, patch_indices_unfilled, clean_indices_unfilled);
    [recon_img, recon_img1, RMSE, zero_pxls] = reconstruct(new_img, img, psz, Y, patch_indices, Y1, patch_indices_unfilled, 'complete');
    
    numel(find(isnan(new_img)));
    numel(find(isnan(recon_img)));

    RMSE
    save(strcat(scene, '_euclidean_RMSE_', num2str(RMSE), '_c', num2str(corrupt_percent(iter)), '_psz', num2str(psz), '_onTheFly', '.mat'), 'recon_img1', 'new_img', 'img', 'RMSE', 'zero_pxls', 'mu', 'lambda', 'l', 'num_patches_each_img', 'A', 'S', 'S1');
end