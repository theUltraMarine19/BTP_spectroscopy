% For German dataset
cd MPRS-161026-1328\
data = fitsread('AspirinANDParacetamol-1_imcmb_oextr.fits', 'image');
data = data + abs(min(data(:)));
data_trans = data';
load('sorting.mat');
cd ..
indices = zeros(1, numel(C));
C(149) = 67;   % HARD-CODE
C(272) = 344;  % HARD-CODE

for i=1:numel(C)
    indices(:, i) = find(C==i);
end

data_final = data_trans(indices(1,:), :);

fileID = fopen('Matrix.txt','r');
formatSpec = '%f';

size_img = 20;
step = 20;
M = data_final;
M_reformed = zeros(size(M,1), floor(size(M,2)/step));
idx = 1;

for j=1:step:size(M,2)-step
    M_reformed(:,idx) = sum(M(:,j:j+step-1), 2);
    idx = idx + 1;
end

ctr = 1;
img = zeros(size_img, size_img, floor(size(M,2)/step));
orig_img = zeros(size_img, size_img, floor(size(M,2)/step));
for i=1:size_img
    for j=1:size_img
        img(j, i, :) = M_reformed(ctr, :);
        orig_img(i, j, :) = M_reformed(ctr, :);
        ctr = ctr + 1;
    end
end

img_avg = img./step;
