function [im] = Qdisp_Dictionary(D)
% function displays the dictionary atoms as blocks.
% each atom would be normalized to [0 255]
% dictionary can be sorted by variance

borderSize = 1;
[m n] = size(D(:,:,1));
k = n;            % number of blocks
bb = sqrt(m);   % blocksize of each rgb_block

num_row = floor(sqrt(k));  % get the number of rows of blocks to be shown 
num_col = floor(sqrt(k));  % get the number of cols of blocks to be shown 
sizeForEachBlk = bb+borderSize;  % get block_size to be shown, usually bb+1
I = zeros(sizeForEachBlk*num_row+borderSize,sizeForEachBlk*num_col+borderSize,3);

% fill all the image with blue color
I(:,:,1) = 0;
I(:,:,2) = 0;
I(:,:,3) = 0;

% fill the elements of dicionary to corresponding block
% first: normailization
for ii = 1:size(D,2)
    D(:,ii,2) = D(:,ii,2) - min(D(:,ii,2));
    if ( max(D(:,ii,2)) )
        D(:,ii,2) = D(:,ii,2)./max(D(:,ii,2));
    end
    D(:,ii,3) = D(:,ii,3) - min(D(:,ii,3));
    if ( max(D(:,ii,3)) )
        D(:,ii,3) = D(:,ii,3)./max(D(:,ii,3));
    end
    D(:,ii,4) = D(:,ii,4) - min(D(:,ii,4));
    if ( max(D(:,ii,4)) )
        D(:,ii,4) = D(:,ii,4)./max(D(:,ii,4));
    end
end

% second: fill dictionary into blocks
counter = 1;
for ii = 1:num_row
    for jj = 1:num_col
        I(borderSize+(ii-1)*sizeForEachBlk+1:ii*sizeForEachBlk,borderSize+(jj-1)*sizeForEachBlk+1:jj*sizeForEachBlk,:)=...
            reshape(D(:,counter,2:4), bb,bb,3);
        counter = counter + 1;
    end
end
im = I;












    
    
    
    
    
