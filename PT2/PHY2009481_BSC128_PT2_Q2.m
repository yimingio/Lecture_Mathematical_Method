%MatLab R2022a
%PHY2009481 Q2
clc;clear;close all;
A=[1,1,0,3;
    2,1,-1,1;
    3,-1,-1,2;
    -1,2,3,-1];
[n,m]=size(A);

[L,U]=func_LU(A);
fprintf("Question 2(a):\n The L matrix is\n");
disp(L);
fprintf("\n The U matrix is\n");
disp(U);
fprintf("-------------------------------------------\nQuestion 2(b):\n");


B=[8;7;14;-7];
Ytest=inv(L)*B;
Xtest=inv(U)*Ytest;


A1=[L,B];
Y(1)=A1(1,n+1)/A1(1,1);
for i=2:n
    sum=0;
    for j=1:i-1
        sum=sum+A1(i,j)*Y(j);
    end
    Y(i)=(A1(i,n+1)-sum)/A1(i,i);
end

Y=Y';

X=[];
A2=[U,Y];
X(n)=A2(n,n+1)/A2(n,n);
for i=n-1:-1:1
    sum=0;
    for j=i+1:n
        sum=sum+A2(i,j)*X(j);
    end
    X(i)=(A2(i,n+1)-sum)/A2(i,i);
end
X=X';
fprintf("(LY=b) After substitution the Y is \n");
disp(Y);
fprintf("(UX=Y) After substitution the X is \n");
disp(X);


function [L,U]=func_LU(A)
[n,m]=size(A);
for i=1:n-1
    %Check and switch to check pivot to be not 0
    p=i;
    while A(p,i)==0 & p <= n-1    
        p=p+1;
    end
    if p==n & A(p,i)==0
        disp("No unique solution found, please check");
        break
    end
    Atemp=A(i,:);
    A(i,:)=A(p,:);
    A(p,:)=Atemp;
    for j=i+1:n
        m(j,i)=A(j,i)/A(i,i);
        A(j,:)=A(j,:)-m(j,i).*A(i,:);
    end
    Ag{i +1}=A;
end





% modify into LU factorization

L2=eye(n);
L1=zeros(4,4);
L1(2:n,1:n-1)=m(2:end,:);
L=L1+L2;
U=A(:,1:n);


end