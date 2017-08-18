% input a and b are quaternion matrices
% output y is also a quaternion matrix, which is the result of the
% multiplication of input a and b
function y = Qmult(a,b)
[ax, ay] = size(a(:,:,1));
[bx, by] = size(b(:,:,1));
y = zeros(ax, by, 4);
y(:,:,1) = a(:,:,1)*b(:,:,1) - a(:,:,2)*b(:,:,2) - a(:,:,3)*b(:,:,3) - a(:,:,4)*b(:,:,4);
y(:,:,2) = a(:,:,1)*b(:,:,2) + a(:,:,2)*b(:,:,1) + a(:,:,3)*b(:,:,4) - a(:,:,4)*b(:,:,3);
y(:,:,3) = a(:,:,1)*b(:,:,3) + a(:,:,3)*b(:,:,1) + a(:,:,4)*b(:,:,2) - a(:,:,2)*b(:,:,4);
y(:,:,4) = a(:,:,1)*b(:,:,4) + a(:,:,4)*b(:,:,1) + a(:,:,2)*b(:,:,3) - a(:,:,3)*b(:,:,2);
