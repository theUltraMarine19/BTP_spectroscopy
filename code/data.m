function X = data(psz, num)
ctr = 1;

X = zeros(psz*psz,1*6*(112-psz+1)*(92-psz+1));   % data matrix whose columns are vectorised yi's
for i = num:num % no. of face folders
    for j = 1:6    % no. of images in each folder
        
        img = imread(strcat('att_faces/s', num2str(i), '/', num2str(j), '.pgm'));
        img = double(img);
        num_patches_each_img = (size(img,1)-psz+1)*(size(img,2)-psz+1);
        patches = zeros(psz*psz, num_patches_each_img);
%         patch_cell = mat2cell(img(:, 3:90), psz*ones(1,size(img,1)/psz), psz*ones(1,(size(img,2)-4)/psz));
        
%         rand_rows = randperm(14);
%         rand_cols = randperm(11);
%         rand_pix = randperm((size(img,1) - psz + 1)*(size(img,2) - psz + 1));
%         rand_idx = rand_pix(1:num_patches_each_img);
        % choose randomly num_patches_each_img patches out of 154 obtained
%         rand_rows = rand_rows(1:num_patches_each_img);
%         rand_cols = rand_cols(1:num_patches_each_img);
        for k = 1:num_patches_each_img
            sz = size(img,2) - psz + 1;
            row = floor((k-1)/sz) + 1;
            col = mod(k + sz -1, sz) + 1;
            tmp_window = img(row:row+psz-1, col:col+psz-1);
            tmp = tmp_window(:);
%             % for each patch, take 80% of the values
%             rand_indices = randperm(64);
%             rand_indices = rand_indices(13:64);
%             patches(:, k) = tmp(rand_indices);
            patches(:, k) = tmp;
            X(:, ctr) = patches(:, k);
            ctr = ctr + 1;
        end
    end
end

        
            
        
        
        
        