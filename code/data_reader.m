size_img = [21 21];                 % dimensions of the image
step = 3;                           % Wavelengths to be skipped
M = csvread('paraffin_skin2.csv');  % CSV file to be read

% Labels for first two columns to be ignored
% spectra = M(1, :);
% spectra(:, 1:2) = []

M(1, :) = [];
M(:, 1:2) = [];
size(M,2)
M_reformed = zeros(size(M,1), floor(size(M,2)/step));
idx = 1;

for j=1:step:size(M,2)-step
    M_reformed(:,idx) = sum(M(:,j:j+step-1), 2);
    idx = idx + 1;
end

ctr = 1;
img = zeros(size_img(1), size_img(2), floor(size(M,2)/step));
orig_img = zeros(size_img(1), size_img(2), size(M,2));
for i=1:size_img(2)
    for j=1:size_img(1)
        img(j, i, :) = M_reformed(ctr, :);
        orig_img(j, i, :) = M(ctr, :);
        ctr = ctr + 1;
    end
end

img_avg = img./step;
