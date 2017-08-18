function [abs_A]=Q_abs(A)
% Absolute value.
P=zeros(size(A(:,:,1)));
P=A(:,:,1).^2+A(:,:,2).^2+A(:,:,3).^2+A(:,:,4).^2;
abs_A=P;
