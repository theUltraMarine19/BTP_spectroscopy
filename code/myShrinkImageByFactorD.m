% Image downsampling code
function [new_img] = myShrinkImageByFactorD(img, d)
% img = imread('circles_concentric.png');
new_img = img(1:d:size(img, 1), 1:d:size(img, 2), :);
daspect([1 1 1]);
imshow(new_img(:,:,1:3))
colorbar
axis on;
end
