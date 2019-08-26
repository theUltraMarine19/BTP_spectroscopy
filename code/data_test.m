function [new_img, X, patch_indices, clean_indices] = data_test(old_image, psz, num_img_patches, corrupt_percent)
new_img = double(old_image);

X = zeros(psz*psz, num_img_patches);
% X = cell(num_img_patches, 1);
% patch_ind_rand = randperm((size(new_img,1) - psz + 1) * (size(new_img,2) - psz + 1));
% patch_indices = patch_ind_rand(1:num_img_patches);
patch_indices = randi([1, (size(new_img,1) - psz + 1) * (size(new_img,2) - psz + 1)], 1, num_img_patches);

% phi = zeros(ceil(psz*psz*0.2), psz*psz);
clean_indices = cell(num_img_patches, 1);
rand_idx = randperm(size(new_img,1)*size(new_img,2));
num_to_corrupt = ceil(corrupt_percent*size(new_img,1)*size(new_img,2));

for i = 1:num_to_corrupt
    row = floor((rand_idx(i)-1)/size(new_img,2)) + 1;
    col = mod(rand_idx(i) + size(new_img,2) - 1, size(new_img,2)) + 1;
    new_img(row, col) = -1.0;  % black out the corrupt pixels (distinguish from original image black pixels)
end

ctr = 1;
sz = size(new_img,2) - psz + 1;
for i = 1:numel(patch_indices)
    row = floor((patch_indices(i)-1)/sz) + 1;
    col = mod(patch_indices(i) + sz -1, sz) + 1;
    img_window = new_img(row:row+psz-1, col:col+psz-1);
    img_vect = img_window(:);
    X(:, ctr) = img_vect;
    clean_indices{ctr} = find(img_vect~=-1.0)';
    ctr = ctr+1;
end
end
