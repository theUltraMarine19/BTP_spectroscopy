function X = data_hyper(img, psz, num_patches_each_img)
ctr = 1;

num_channels = size(img, 3);
X = zeros(psz*psz*num_channels, num_patches_each_img);   % data matrix whose columns are vectorised yi's
        
sz_horiz = size(img,2);

rand_pix = randperm((size(img,1) - psz + 1)*(sz_horiz - psz + 1));
rand_idx = rand_pix(1:num_patches_each_img);
for k = 1:num_patches_each_img
    sz = sz_horiz - psz + 1;
    row = floor((rand_idx(k)-1)/sz) + 1;
    col = mod(rand_idx(k) + sz -1, sz) + 1;
    tmp_window = img(row:row+psz-1, col:col+psz-1, :);
    tmp = tmp_window(:);
    X(:, ctr) = tmp;
    ctr = ctr + 1;
end

        
            
        
        
        
        