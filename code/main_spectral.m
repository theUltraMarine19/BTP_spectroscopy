function [] = main_spectral(img, scene, D2)
lambda = 2e-4;
mu = 3e-5;
num_patches_each_img = [18*18, 5000, 8000, 10500];
corrupt_percent = [0, 0.2 0.5 0.8];
psz = 4;
l = 12;

%% Adding noise here
% img1 = poissrnd(img, 'poisson');

%% Main Script

for iter = 1:1
    for refine=1:3
    [new_img, X, X1, X_super, patch_indices, patch_indices_unfilled, clean_indices, clean_indices_unfilled, super_indices] = data_test_hyper(img, psz, num_patches_each_img(iter), corrupt_percent(iter), 'complete', 'random', 'DC');
    [Y, A, S] = spectral_sep(X, lambda, l, psz, mu, num_patches_each_img(iter), patch_indices, clean_indices, 2, D2);
%     [Y1, S1] = test_fast(A, X1, lambda, patch_indices_unfilled, clean_indices_unfilled);
    [recon_img, recon_img1, RMSE, zero_pxls] = reconstruct(new_img, img, psz, Y, patch_indices, [], [], 'partial');
    
    img = recon_img1;

    save(strcat(scene, '_spectral', '_psz', num2str(psz), '_onTheFly', '.mat'), 'recon_img1', 'new_img', 'img', 'RMSE', 'mu', 'lambda', 'l', 'num_patches_each_img', 'A', 'S');
    end
end