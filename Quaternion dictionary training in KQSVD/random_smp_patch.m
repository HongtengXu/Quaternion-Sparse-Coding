function [X] = random_smp_patch(IMG_PATH, type, patch_size, num_patch)
% randomly sample image patches
img_dir = dir(fullfile(IMG_PATH, type));

X = [];

img_num = length(img_dir);
nper_img = zeros(1,img_num);   % compute how many pixels in per image 

for ii = 1:img_num
    im = imread(fullfile(IMG_PATH, img_dir(ii).name));
    nper_img(ii) = prod(size(im(:,:,1)));
end

nper_img = floor(nper_img*num_patch/sum(nper_img));% compute num of patches each pic

for ii = 1:img_num
    tmp_num = nper_img(ii);
    im = im2double(imread(fullfile(IMG_PATH, img_dir(ii).name)));  % change into double type
    [tmp_X] = sample_patches(im, patch_size, tmp_num);
    X = [X, tmp_X];
end

% patch_path = ['Training result/' num2str(num_patch) '_patches' '.mat'];
% save(patch_path, 'X');

    