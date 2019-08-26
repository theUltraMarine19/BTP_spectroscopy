function [new_img, new_img_interp, X, X1, X_interp, patch_indices, patch_indices_unfilled, clean_indices, clean_indices_unfilled] = data_gmm(old_image, psz, num_img_patches, corrupt_percent, corrupt_percent_spectral, sampling, corruption_pattern)
% corrupt_percent : spatial undersampling
% corrupt_percent_spectral : spectral undersampling (Set it to 1 for full blackout
% sampling : random / complete
% corruption_pattern : random / (periodic for super-res)

new_img = double(old_image);
mask = zeros(size(old_image));
tracker = zeros(size(old_image,1), size(old_image,2), size(old_image, 3));

num_channels = size(new_img, 3);
X = zeros(psz*psz*psz, num_img_patches);
sz_horiz = size(new_img,2);

if strcmp(sampling, 'random') == 1  % select patches randomly
    patch_indices = randi([1, (size(new_img,1) - psz + 1) * (sz_horiz - psz + 1) * (num_channels - psz + 1)], 1, num_img_patches);
else                                % raster scan for patch selection
    if num_img_patches > (size(new_img,1) - psz + 1) * (sz_horiz - psz + 1) * (num_channels - psz+1)
        return;
    end
    patch_indices = 1:num_img_patches;
end

clean_indices = cell(num_img_patches, 1);
rand_idx = randperm(size(new_img,1)*sz_horiz);
num_to_corrupt = ceil(corrupt_percent*size(new_img,1)*sz_horiz);
num_to_corrupt_spectral = ceil(corrupt_percent_spectral*num_channels);

if strcmp(corruption_pattern, 'random') == 1
    for i = 1:num_to_corrupt
        row = floor((rand_idx(i)-1)/sz_horiz) + 1;
        col = mod(rand_idx(i) + sz_horiz - 1, sz_horiz) + 1;
        idx_spectral = randperm(num_channels, num_to_corrupt_spectral);
        new_img(row, col, idx_spectral) = -1.0;  % black out the corrupt pixels (distinguish from original image black pixels)
        mask(row, col, idx_spectral) = 1;
    end
% Structural not yet done for GMMs
else % periodic structural corruption for super-res 
    new_img = zeros(size(old_image)) - 1.0;
    new_img(1:corrupt_percent:size(old_image,1), 1:corrupt_percent:size(old_image,2), 1:corrupt_percent_spectral:num_channels) = double(old_image(1:corrupt_percent:size(old_image,1), 1:corrupt_percent:size(old_image,2), :));
end

% ctr = 1;
dims = [size(new_img,1) - psz + 1, size(new_img,2) - psz + 1, num_channels-psz+1];
% sz = sz_horiz - psz + 1;
for i = 1:numel(patch_indices)
%     row = floor((patch_indices(i)-1)/sz) + 1;
%     col = mod(patch_indices(i) + sz -1, sz) + 1;
    [row, col, sp_idx] = ind2sub(dims, patch_indices(i));
    
    img_window = new_img(row:row+psz-1, col:col+psz-1, sp_idx:sp_idx+psz-1);
    tracker(row:row+psz-1, col:col+psz-1, sp_idx:sp_idx+psz-1) = ones(psz, psz, psz);
    img_vect = img_window(:);
    X(:, patch_indices(i)) = img_vect;
    clean_indices{patch_indices(i)} = find(img_vect~=-1.0)'; % clean indices in each psz X psz X psz patch
%     ctr = ctr+1;

end

% X1 is the matrix consisting of those patches which were not sampled (for hole filling) %
tracker(:,(size(new_img,2) - psz + 1)) = 0;     % Right border
tracker((size(new_img,1) - psz + 1), :) = 0;    % Bottom border
unfilled = find(tracker(1:(size(new_img,1) - psz + 1),1:(size(new_img,2) - psz + 1),1:num_channels-psz+1)==0);
X1 = zeros(psz*psz*psz, numel(unfilled));
clean_indices_unfilled = cell(numel(unfilled), 1);

ctr1 = 1;
patch_indices_unfilled = zeros(1, numel(unfilled));

for i=1:numel(unfilled)
    unfilled(i);
%     col = floor((unfilled(i)-1)/(size(new_img,1)-psz+1)) + 1;
%     row = mod(unfilled(i) + size(new_img,1) - psz, (size(new_img,1)-psz+1)) + 1;
    [row, col, sp_idx] = ind2sub(dims, unfilled(i));
%     patch_indices_unfilled(:,i) = (row-1)*(sz_horiz-psz+1) + col;
    patch_indices_unfilled(:,i) = sub2ind(dims, row, col, sp_idx);
    img_window = new_img(row:row+psz-1, col:col+psz-1, sp_idx:sp_idx+psz-1);
    img_vect = img_window(:);
    X1(:, ctr1) = img_vect;
    clean_indices_unfilled{ctr1} = find(img_vect~=-1.0)';
    ctr1 = ctr1 + 1;
end

new_img_interp = zeros(size(new_img,1), size(new_img,2), size(new_img,3));
for ch=1:size(new_img,3)
    new_img_interp(:,:,ch) = regionfill(new_img(:,:,ch), mask(:,:,ch));  
end

X_interp = zeros(psz*psz*psz, (size(new_img,1)-psz+1)*(size(new_img,2)-psz+1)*(num_channels-psz+1));
for i = 1:numel(patch_indices)
    [row, col, sp_idx] = ind2sub(dims, patch_indices(i));
    img_window = new_img_interp(row:row+psz-1, col:col+psz-1, sp_idx:sp_idx+psz-1);
    img_vect = img_window(:);
    X_interp(:, patch_indices(i)) = img_vect;
end