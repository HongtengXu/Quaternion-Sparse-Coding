function X = Qpinv(A,varargin)
%PINV   Pseudoinverse.
%   X = PINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS(class(A)).
%
%   PINV(A,TOL) uses the tolerance TOL instead of the default.
%
%   Class support for input A: 
%      float: double, single
%
%   See also RANK.



[m,n,dim] = size(A);

if n > m
   %X = pinv(A',varargin{:})';
   B(:,:,1)=A(:,:,1)';
   B(:,:,2)=A(:,:,2)';
   B(:,:,3)=A(:,:,3)';
   B(:,:,4)=A(:,:,4)';
   XX=Qpinv(B,varargin{:});
   X(:,:,1)=XX(:,:,1)';
   X(:,:,2)=XX(:,:,2)';
   X(:,:,3)=XX(:,:,3)';
   X(:,:,4)=XX(:,:,4)';
else
   [Uq,Sq,Vq] = qsvd(A);
   if m > 1, s = diag(Sq);
      elseif m == 1, s = Sq(1);
      else s = 0;
   end
   if nargin == 2
      tol = varargin{1};
   else
      tol = max(m,n) * eps(max(s));
   end
   r = sum(s > tol);
   if (r == 0)
      X = [];
   else
      s = diag(ones(r,1)./s(1:r));
      %X = V(:,1:r)*s*U(:,1:r)';
      [aa,bb] = size(s);
      ss = zeros(aa,bb,4);
      ss(:,:,1) = s;
      X1 = Qmult(Vq(:,1:r,:), ss);
      Uq1 = Uq(:,1:r,:);
      Uq2(:,:,1) = Uq1(:,:,1)';
      Uq2(:,:,2) = -Uq1(:,:,2)';
      Uq2(:,:,3) = -Uq1(:,:,3)';
      Uq2(:,:,4) = -Uq1(:,:,4)';
      X = Qmult(X1, Uq2);
   end
end
