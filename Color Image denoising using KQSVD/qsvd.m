function [Uq,Sq,Vq] = qsvd(Qq)
%
% GENERAL DESCRIPTION
% Quaternionic Singular Value Decomposition
% [Uq,Sq,Vq] = qsvd(Qq) produces a diagonal matrix Sq of the same dimension as the real 
% part of a quaternionic matrix Qq, with nonnegative diagonal elements in decreasing 
% order, and quaternionic unitary matrices Uq and Vq so that Qq = Uq*Sq*Vq'.
%
% ABOUT
% Create:               14 Mar. 2008
% Last Update:       14 Mar. 2008
% Reversion:           1.0
%
% AUTHOR            Lu, Wei
%
% AFFILIATION	Institute of Image Communication and Information Processing,
%                           Shanghai Jiao Tong University, China
% 
% Mailto:               learza2008@gmail.com 
% 
% REFERENCE       F. Zhang, "Quaternions and matrices of Quaternions," in
%                           Linear Algebra and Its Applications, 1997, pp. 21-57.

% initialization
%sz = size(Qq(:,:,1));
%Uq = zeros(sz(1),sz(1),4);
%Vq = zeros(sz(2),sz(2),4);

% computer svd of the equivalent complex matrix of Qq
Ac = Qq(:,:,1)+i*Qq(:,:,2);
Bc = Qq(:,:,3)+i*Qq(:,:,4);
Cc = [Ac Bc;-conj(Bc) conj(Ac)];
[Uc S Vc] = svd(Cc, 'econ');

% use the relationship between Cc and Qq to obtain Uq, Sq, and Vq
Uq(:,:,1) = real(Uc(1:end/2,1:2:end));
Uq(:,:,2) = imag(Uc(1:end/2,1:2:end));
Uq(:,:,3) = real(-conj(Uc(end/2+1:end,1:2:end)));
Uq(:,:,4) = imag(-conj(Uc(end/2+1:end,1:2:end)));

Vq(:,:,1) = real(Vc(1:end/2,1:2:end));
Vq(:,:,2) = imag(Vc(1:end/2,1:2:end));
Vq(:,:,3) = real(-conj(Vc(end/2+1:end,1:2:end)));
Vq(:,:,4) = imag(-conj(Vc(end/2+1:end,1:2:end)));


[mm,nn]=size(S);
ii = 1;
j = 1;
Sq = zeros(min(mm,nn)/2,min(mm,nn)/2);
while ii <= min(mm,nn)
Sq(j,j) = S(ii,ii);
j = j+1;
ii = ii+2;
end;



