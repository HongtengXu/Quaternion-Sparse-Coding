%==========================================================================
% Authors:Yi Xu, Licheng Yu, Hongteng Xu, Hao Zhang and Truong Nguyen  
% TItle:Vector Sparse Representation of Color Image Using Quaternion Matrix Analysis
% IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 24, NO. 4, PP.1315-1329, APRIL 2015,
% Codes for color image reconstruction using K-QSVD

clear all;clc;close all;
IMG_PATH = 'Dataset classifying/animal/';

Im = im2double(imread([IMG_PATH 'frog.ppm']));
Reduce_DC = 1;
[m n] = size(Im(:,:,1)); % size of input image
bb = 8;                  % size of each block
step = 8;                % sliding distance
sparsity = 5;

[blocks,idx] = Q_im2col(Im, bb, step);

load('Training result/dict_256_atoms.mat');
D = Dictionary;
tt = cputime;
for jj=1:10000:size(blocks,2)
    jump_size = min(jj+10000-1,size(blocks,2));
    if (Reduce_DC)
        mean_R = repmat(mean(blocks(:,jj:jump_size,2)),size(blocks(:,jj:jump_size,:),1),1);
        mean_G = repmat(mean(blocks(:,jj:jump_size,3)),size(blocks(:,jj:jump_size,:),1),1);
        mean_B = repmat(mean(blocks(:,jj:jump_size,4)),size(blocks(:,jj:jump_size,:),1),1);
        blocks(:,jj:jump_size,2) = blocks(:,jj:jump_size,2) - mean_R;
        blocks(:,jj:jump_size,3) = blocks(:,jj:jump_size,3) - mean_G;
        blocks(:,jj:jump_size,4) = blocks(:,jj:jump_size,4) - mean_B;
    end    
    
    Coefs = QOMP(D, blocks(:,jj:jump_size,:),sparsity);
    
    if (Reduce_DC)
        blocks(:,jj:jump_size,:) = Qmult(D,Coefs);
        blocks(:,jj:jump_size,2) = blocks(:,jj:jump_size,2) + mean_R;
        blocks(:,jj:jump_size,3) = blocks(:,jj:jump_size,3) + mean_G;
        blocks(:,jj:jump_size,4) = blocks(:,jj:jump_size,4) + mean_B;
    else
        blocks(:,jj:jump_size) = Qmult(D,Coefs);
    end
end

Reconst_IM = zeros(size(Im));
Weight = zeros(size(Im));

for ii=1:3
    tmp_Blks = blocks(:,:,1+ii);
    [rows,cols] = ind2sub(size(Reconst_IM(:,:,ii))-bb+1,idx); % revise required
    count = 1;
    for jj = 1:length(rows)
        row = rows(jj);
        col = cols(jj);
        blk = reshape(tmp_Blks(:,count,:),[bb bb]);
        Reconst_IM(row:row+bb-1,col:col+bb-1,ii) = Reconst_IM(row:row+bb-1,col:col+bb-1,ii)+blk;
        Weight(row:row+bb-1,col:col+bb-1,ii) = Weight(row:row+bb-1,col:col+bb-1,ii)+1;
        count = count+1;
    end
end
Reconst_IM = Reconst_IM./Weight;

tt = cputime-tt;
fprintf('using time: %ds\n',tt)

PSNR = 20*log10(1/sqrt(mean((Reconst_IM(:)-Im(:)).^2)));
imshow(Reconst_IM,[]);



    





