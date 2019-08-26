function [] = main_fly_poisson(img, scene)
lambda = 1e-4;
mu = 1.5e-5;
num_patches_each_img = [4200 4800 5400 6000];
corrupt_percent = [0 0.2 0.5 0.8];
psz = 8;
l = 20;

%% Adding noise here
img1 = poissrnd(img./100.0);

%% Main Script

for iter = 1:1
    [new_img, X, X1, patch_indices, patch_indices_unfilled, clean_indices, clean_indices_unfilled] = data_test_hyper(img1, psz, num_patches_each_img(iter), corrupt_percent(iter));
    [Y, A, S] = compute_fly_batch_poisson(X, lambda, l, psz, mu, num_patches_each_img(iter), patch_indices, clean_indices, 18);
    [Y1, S1] = test_fast_poisson(A, X1, lambda, patch_indices_unfilled, clean_indices_unfilled);
    [recon_img, recon_img1, RMSE, zero_pxls] = reconstruct(new_img, img, psz, Y, patch_indices, Y1, patch_indices_unfilled);
    
    numel(find(isnan(new_img)));
    numel(find(isnan(recon_img)));

    RMSE
    save(strcat(scene, '_poisson_RMSE_', num2str(RMSE), '_c', num2str(corrupt_percent(iter)), '_psz', num2str(psz), '_onTheFly', '.mat'), 'recon_img1', 'new_img', 'img', 'RMSE', 'zero_pxls', 'mu', 'lambda', 'l', 'num_patches_each_img');
end