% Generalized reconstruction
function [recon_img, recon_img1, RMSE, zero_pxls] = reconstruct(test_img, old_img, psz, Y, patch_indices, Y1, patch_indices_unfilled, mode)
recon_img = double(test_img);
recon_img1 = double(zeros(size(test_img)));
old_img = double(old_img);
counts1 = zeros(size(test_img));
counts = ones(size(test_img));
mask = zeros(size(test_img));

num_channels = size(test_img,3);

% Reconstruct the image
for i = 1:numel(patch_indices)
    sz = size(test_img,2) - psz + 1;
    row = floor((patch_indices(i)-1)/sz) + 1;
    col = mod(patch_indices(i) + sz -1, sz) + 1;
    tmp = reshape(Y(:, i), [psz psz num_channels]);
    recon_img(row:row+psz-1, col:col+psz-1, :) = recon_img(row:row+psz-1, col:col+psz-1, :) + tmp;
    recon_img1(row:row+psz-1, col:col+psz-1, :) = recon_img1(row:row+psz-1, col:col+psz-1, :) + tmp;
    counts(row:row+psz-1, col:col+psz-1, :) = counts(row:row+psz-1, col:col+psz-1, :) + ones(psz, psz, num_channels);
    counts1(row:row+psz-1, col:col+psz-1, :) = counts1(row:row+psz-1, col:col+psz-1, :) + 1;
    mask(row:row+psz-1, col:col+psz-1, :) = mask(row:row+psz-1, col:col+psz-1, :) + 1;
end

if strcmp(mode, 'complete') == 1
    for i = 1:numel(patch_indices_unfilled)
        sz = size(test_img,2) - psz + 1;
        row = floor((patch_indices_unfilled(i)-1)/sz) + 1;
        col = mod(patch_indices_unfilled(i) + sz -1, sz) + 1;
        tmp = reshape(Y1(:, i), [psz psz num_channels]);
        recon_img(row:row+psz-1, col:col+psz-1, :) = recon_img(row:row+psz-1, col:col+psz-1, :) + tmp;
        recon_img1(row:row+psz-1, col:col+psz-1, :) = recon_img1(row:row+psz-1, col:col+psz-1, :) + tmp;
        counts(row:row+psz-1, col:col+psz-1, :) = counts(row:row+psz-1, col:col+psz-1, :) + ones(psz, psz, num_channels);
        counts1(row:row+psz-1, col:col+psz-1, :) = counts1(row:row+psz-1, col:col+psz-1, :) + 1;
        mask(row:row+psz-1, col:col+psz-1, :) = mask(row:row+psz-1, col:col+psz-1, :) + 1;
    end
end

% sum_error = zeros(64,1);
% sum_error = 0;
zero_pxls = numel(find(counts1==0));

recon_img = recon_img./counts;
recon_img1 = recon_img1./counts1;
numel(find(isnan(recon_img)));

% for i = 1:numel(patch_indices)
%     sz = size(test_img,2) - psz + 1;
%     row = floor((patch_indices(i)-1)/sz) + 1;
%     col = mod(patch_indices(i) + sz -1, sz) + 1;
%     diff = recon_img(row:row+psz-1, col:col+psz-1, :) - old_img(row:row+psz-1, col:col+psz-1, :);
%     true = old_img(row:row+psz-1, col:col+psz-1, :);
%     sum_error = sum_error + norm(diff(:), 2)./norm(true(:), 2);
% end

% sum_error = sum_error/numel(patch_indices);

mask(mask~=0) = 1;
recon_img1(isnan(recon_img1)) = 0;

RMSE = norm(mask(:).*(recon_img1(:) - old_img(:)),2)/norm(mask(:).*old_img(:),2);