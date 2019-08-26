% operations for ORL face database : Dictionary learning, now utilized for GMM
function [recon_img1, mask] = main(img2, fname)
for num=2:5
    fname = strcat('s', num2str(num), '_8');
    img2 = imread(strcat('att_faces\s', num2str(num), '\8.pgm'));
    lambda_train = 2e-4;
    lambda_test = 4e-4;
    mu = 5e-5;
    num_patches_each_img = (112-8+1)*(92-8+1);
    corrupt_percent = [0.2 0.5 0.8];
    psz = 8;
    k = 25;
    
    
    for iter=1:3
        X = data(psz,num);
        % [A, S] = compute(X, lambda_train, 35, mu, psz, num_patches_each_img);
        [new_img, X_test, patch_indices, clean_indices] = data_test(img2, psz, num_patches_each_img, corrupt_percent(iter));
        avg_intensity = sum(new_img(:))/numel(new_img);
        noise_var = (0.01 * avg_intensity)^2;
        % [Y,S] = test_fast(A, X_test, lambda_test, patch_indices, clean_indices);
        [label, model, llh] = mixGaussEm(X, k);
        Z = gmm_estimate(X_test, model, psz*psz, patch_indices, clean_indices, noise_var)
        [recon_img, recon_img1, RMSE, mask] = reconstruct(new_img, img2, psz, Z, patch_indices, [], [], 'complete');
        % recon_img = reconstruct(old_img, old_img, 8, Y, patch_indices) % WON'T WORK because the Y will be for corrupted image only

        % numel(find(isnan(new_img)));
        % numel(find(isnan(recon_img)));
        % 
        % maxc = max(recon_img1(:));
        % minc = min(recon_img1(:));
        % maxi = max(new_img(:));
        % mini = min(new_img(:));
        % maxo = max(img2(:));
        % mino = min(img2(:));

        % for i = 1:floor(size(img1,3)/3)
        %     img1 = mat2gray(img1);
        %     new_img = mat2gray(new_img);
        %     recon_img1 = mat2gray(recon_img1);

        %     figure;
        %     subplot(1,3,1);
        %     imshow(img2(:,:,3*i-2:3*i), [mino, maxo]);
            % imwrite(mat2gray(recon_img1), 's23_orig.png');
        %     subplot(1,3,2);
        %     imshow(new_img(:,:,3*i-2:3*i), [mini, maxi]);
            % imwrite(mat2gray(recon_img1), 's23_noisy.png');
            % subplot(2,2,3);
            % imshow(mat2gray(recon_img));
        %     subplot(1,3,3);
        %     imshow(recon_img1(:,:,3*i-2:3*i), [minc, maxc]);
            imwrite(mat2gray(recon_img1), strcat(fname, '_', 'c', num2str(corrupt_percent(iter)), '_psz', num2str(psz), '_RMSE', num2str(RMSE), '.png'));
    
    end
    % RMSE
end