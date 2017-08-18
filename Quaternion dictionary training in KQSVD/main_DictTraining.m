%==========================================================================
% Authors:Yi Xu, Licheng Yu, Hongteng Xu, Hao Zhang and Truong Nguyen  
% TItle:Vector Sparse Representation of Color Image Using Quaternion Matrix Analysis
% IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 24, NO. 4, PP.1315-1329, APRIL 2015,
% Codes for dictionary training by K-QSVD


clear all;clc;close all;
IMG_PATH = 'Dataset classifying/flower2/';

dict_size = 126;       % dictionary size
L = 5;                 % number of non-zero coeficients for each representation
Out_iter = 10;          % number of iterations of KSVD algorithm
patch_size = 8;        % image patch size
num_patch = 10002;     % number of patches

% randomly sample image patches
[temp_X] = random_smp_patch(IMG_PATH, '*.bmp', patch_size, num_patch); % temp_X is RGB cascaded in column 

hh = size(temp_X,1)/3;
X = zeros( hh, size(temp_X,2), 4 ); % sampled data in Quaternion Matrix Form
X(:,:,2) = temp_X( 1:hh,:); % X(:,:,2) = R channel of sampled data
X(:,:,3) = temp_X( hh+1:hh*2,:); % X(:,:,3) = G channel of sampled data
X(:,:,4) = temp_X( hh*2+1:hh*3,:); % X(:,:,4) = B channel of sampled data

% train the dictionary using KSVD
Reduce_DC = 1;
param.K = dict_size;
param.L = L;
param.numIteration = Out_iter;
param.errorFlag = 0;
param.preserveDCAtom = 0;
param.InitializationMethod = 'DataElements';
param.displayProgress = 1;
if (Reduce_DC)
    mean_R = repmat(mean(X(:,:,2)),size(X,1),1);
    mean_G = repmat(mean(X(:,:,3)),size(X,1),1);
    mean_B = repmat(mean(X(:,:,4)),size(X,1),1);
    X(:,:,2) = X(:,:,2) - mean_R;
    X(:,:,3) = X(:,:,3) - mean_G;
    X(:,:,4) = X(:,:,4) - mean_B;
end

[Dictionary, output] = K_QSVD(X,param);

y = 1:10;
plot(y,output.totalerr2)

figure
[Im_dict] = Qdisp_Dictionary(Dictionary); 
imshow(Im_dict,[]);
%imwrite(Im_dict,'Dic.png');
dict_path = ['Training result/' 'dict_' num2str(dict_size) '_atoms' '.mat'];
save(dict_path, 'Dictionary');



