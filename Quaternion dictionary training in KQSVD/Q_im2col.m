function [blocks,idx] = Q_im2col(I, bb, step)
% Transform input image blocks to column vector
idxMat = zeros(size(I(1:end-bb+1,1:end-bb+1,1)));
idxMat([1:step:end-1,end], [1:step:end-1,end]) = 1;
idx = find(idxMat);

[rows, cols] = ind2sub(size(idxMat), idx);
blocks = zeros(bb*bb, length(idx), 4);

for ii=1:length(rows)
    row = rows(ii);
    col = cols(ii); 
    blk_r = I(row:row+bb-1,col:col+bb-1,1);
    blocks(:,ii,2) = blk_r(:);
    blk_g = I(row:row+bb-1,col:col+bb-1,2);
    blocks(:,ii,3) = blk_g(:);
    blk_b = I(row:row+bb-1,col:col+bb-1,3);
    blocks(:,ii,4) = blk_b(:);
end
    



