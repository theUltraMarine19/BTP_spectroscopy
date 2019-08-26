% Used for patch assembly on data, generating undersampled images, etc.
function [new_img, X, X1, X_super, patch_indices, patch_indices_unfilled, clean_indices, clean_indices_unfilled, super_patch_indices] = data_test_hyper(old_image, psz, num_img_patches, corrupt_percent, sampling, corruption_pattern, init_img_pattern)
new_img = double(old_image);
tracker = zeros(size(old_image,1), size(old_image,2));

num_channels = size(new_img, 3);
X = zeros(psz*psz*num_channels, num_img_patches);
sz_horiz = size(new_img,2);

if strcmp(sampling, 'random') == 1
    patch_indices = randi([1, (size(new_img,1) - psz + 1) * (sz_horiz - psz + 1)], 1, num_img_patches);
else
    if num_img_patches > (size(new_img,1) - psz + 1) * (sz_horiz - psz + 1)
        return;
    end
    patch_indices = 1:num_img_patches;
end

clean_indices = cell(num_img_patches, 1);
rand_idx = randperm(size(new_img,1)*sz_horiz);
num_to_corrupt = ceil(corrupt_percent*size(new_img,1)*sz_horiz);

if strcmp(corruption_pattern, 'random') == 1
    for i = 1:num_to_corrupt
        row = floor((rand_idx(i)-1)/sz_horiz) + 1;
        col = mod(rand_idx(i) + sz_horiz - 1, sz_horiz) + 1;
        new_img(row, col, :) = -1.0;  % black out the corrupt pixels (distinguish from original image black pixels)
    end
else
    new_img = zeros(size(old_image)) - 1.0;
    new_img(1:corrupt_percent:size(old_image,1), 1:corrupt_percent:size(old_image,2), :) = double(old_image(1:corrupt_percent:size(old_image,1), 1:corrupt_percent:size(old_image,2), :));
end
ctr = 1;
sz = sz_horiz - psz + 1;
for i = 1:numel(patch_indices)
    row = floor((patch_indices(i)-1)/sz) + 1;
    col = mod(patch_indices(i) + sz -1, sz) + 1;
    
    img_window = new_img(row:row+psz-1, col:col+psz-1, :);
    tracker(row:row+psz-1, col:col+psz-1) = ones(psz, psz);
    img_vect = img_window(:);
    X(:, ctr) = img_vect;
    
    clean_indices{ctr} = find(img_vect~=-1.0)';
    ctr = ctr+1;
end

% Get a smaller portion of the image (sampled pixels) for super-res and corr. matrix
if strcmp(init_img_pattern, 'super') == 1
    new_img_super = new_img(1:corrupt_percent:size(old_image,1), 1:corrupt_percent:size(old_image,2), :);
    X_super = zeros(psz*psz*num_channels, (size(new_img_super,1)-psz+1)*(size(new_img_super,2)-psz+1));
    sz_super = size(new_img_super,2)-psz+1;
    for i = 1:(size(new_img_super,1)-psz+1)*(size(new_img_super,2)-psz+1)
        row = floor((i-1)/sz_super) + 1;
        col = mod(i + sz_super -1, sz_super) + 1;
        img_window = new_img_super(row:row+psz-1, col:col+psz-1, :);
        img_vect = img_window(:);
        X_super(:, i) = img_vect;
    end
    super_patch_indices = 1:1:(size(new_img_super,1)-psz+1)*(size(new_img_super,2)-psz+1);
    
% For interpolation, interpolate sub-sampled image and form corr. matrix    
elseif strcmp(init_img_pattern, 'interp') == 1
    new_img_super = new_img(1:corrupt_percent:size(old_image,1), 1:corrupt_percent:size(old_image,2), :);
    new_img_interp = zeros(size(new_img,1), size(new_img,2), size(new_img,3));
    for ch=1:size(new_img,3)
        new_img_interp(:,:,ch) = imresize(new_img_super(:,:,ch), [size(new_img,1), size(new_img,2)]);  
    end
    X_super = zeros(psz*psz*num_channels, (size(new_img,1)-psz+1)*(size(new_img,2)-psz+1));
    sz_super = size(new_img,2)-psz+1;
    for i = 1:(size(new_img_interp,1)-psz+1)*(size(new_img_interp,2)-psz+1)
        row = floor((i-1)/sz_super) + 1;
        col = mod(i + sz_super -1, sz_super) + 1;
        img_window = new_img_interp(row:row+psz-1, col:col+psz-1, :);
        img_vect = img_window(:);
        X_super(:, i) = img_vect;
    end
    super_patch_indices = 1:1:(size(new_img_interp,1)-psz+1)*(size(new_img_interp,2)-psz+1);
    
% Above methods are only super-res heuristics
else
    X_super = [];
    super_patch_indices = [];
end
    

% X1 is the matrix consisting of those patches which were not sampled (for hole filling) %
tracker(:,(size(new_img,2) - psz + 1)) = 0;
tracker((size(new_img,1) - psz + 1), :) = 0;
unfilled = find(tracker(1:(size(new_img,1) - psz + 1),1:(size(new_img,2) - psz + 1))==0);
X1 = zeros(psz*psz*num_channels, numel(unfilled));
clean_indices_unfilled = cell(numel(unfilled), 1);

ctr1 = 1;
patch_indices_unfilled = zeros(1, numel(unfilled));

for i=1:numel(unfilled)
    unfilled(i);
    col = floor((unfilled(i)-1)/(size(new_img,1)-psz+1)) + 1;
    row = mod(unfilled(i) + size(new_img,1) - psz, (size(new_img,1)-psz+1)) + 1;
    patch_indices_unfilled(:,i) = (row-1)*(sz_horiz-psz+1) + col;
    img_window = new_img(row:row+psz-1, col:col+psz-1, :);
    img_vect = img_window(:);
    X1(:, ctr1) = img_vect;
    clean_indices_unfilled{ctr1} = find(img_vect~=-1.0)';
    ctr1 = ctr1 + 1;
end