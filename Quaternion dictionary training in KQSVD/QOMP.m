function [A]=QOMP(D,X,L) 
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: 
%       D - the dictionary (its columns MUST be normalized).
%       X - the signals to represent
%       L - the max. number of coefficients for each signal.
% output arguments: 
%       A - sparse coefficient matrix.
%=============================================
[n,P,dim]=size(X);
[n,K,dim]=size(D);
A = zeros(K,P,4);   
    DT(:,:,1)=D(:,:,1)';
    DT(:,:,2)=-D(:,:,2)';
    DT(:,:,3)=-D(:,:,3)';
    DT(:,:,4)=-D(:,:,4)';
for k=1:1:P,
    a=[];
    x=X(:,k,:);
    residual=x;
    indx=zeros(L,1);

    for j=1:1:L
        proj=Qmult(DT,residual);
        [maxVal,pos]=max(Q_abs(proj));%find maximum correlation index
        pos=pos(1);
        indx(j)=pos;
        a = Qmult(Qpinv(D(:,indx(1:j),:)),x);
        residual=x-Qmult(D(:,indx(1:j),:),a);
        if sum(residual(:,:,1).^2+residual(:,:,2).^2+residual(:,:,3).^2+residual(:,:,4).^2) < 1e-9
            break;
        end
    end;
    temp=zeros(K,1,4);
    temp(indx(1:j),1,:)=a;
    A(:,k,:)=temp;
end;
return;
