function [A]=QOMPerr(D,X,errorGoal) 
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: 
%       D - the dictionary (its columns MUST be normalized).
%       X - the signals to represent
%       errorGoal - the maximal allowed representation error for
%                  each siganl.
% output arguments: 
%       A - sparse coefficient matrix.
%=============================================
[n,P,dim]=size(X);
[n,K,dim]=size(D);
E2 = errorGoal^2*2.8*n;%%%%%%%%%%%%
maxNumCoef = n/4;
A = zeros(K,P,4);
    DT(:,:,1)=D(:,:,1)';
    DT(:,:,2)=-D(:,:,2)';
    DT(:,:,3)=-D(:,:,3)';
    DT(:,:,4)=-D(:,:,4)';
    D_R=D(:,:,2);D_G=D(:,:,3);D_B=D(:,:,4);D_O=D(:,:,1);
    Q_ones=ones(n,n);
for k=1:1:P,
    a=[];
    x=X(:,k,:);
    residual=x;
    indx=[];%%%%%%%%%%%%%
    curResNorm2 = sum(residual(:,:,1).^2+residual(:,:,2).^2+residual(:,:,3).^2+residual(:,:,4).^2);%%%%%%%%%%%
    j = 0;
    while curResNorm2>E2 & j<maxNumCoef,
        j = j + 1;
        proj=Qmult(DT,residual);
%         E_mult=D_R'*Q_ones*residual(:,:,2)+D_G'*Q_ones*residual(:,:,3)+D_B'*Q_ones*residual(:,:,4);
        E_mult=D_R'*Q_ones*residual(:,:,2)+D_G'*Q_ones*residual(:,:,3)+D_B'*Q_ones*residual(:,:,4)+D_O'*Q_ones*residual(:,:,1);
        new_proj=Q_abs(proj)+(0.0/n)*E_mult;
        [maxVal,pos]=max(new_proj);%find maximum correlation index   
        pos=pos(1);
        indx(j)=pos;
        a = Qmult(Qpinv(D(:,indx(1:j),:)),x);
        residual=x-Qmult(D(:,indx(1:j),:),a);
        curResNorm2 = sum(residual(:,:,1).^2+residual(:,:,2).^2+residual(:,:,3).^2+residual(:,:,4).^2);%%%%%%%%%%%
    end;
    temp=zeros(K,1,4);
    if(isempty(indx)~=1)
    temp(indx(1:j),1,:)=a;  
    end;
    A(:,k,:)=temp;

end;
return;
