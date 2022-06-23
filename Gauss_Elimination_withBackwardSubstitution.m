%Gauss Elimination with backward substitution
%matlab

clc;clear;close all;
A=[1,-1,2,-1;
    2,-2,3,-3;
    1,1,1,0;
    1,-1,4,3];
B=[-8;
    -20;
    -2;
    4];
Aug=[A B];
N=length(B);
X=zeros(N,1);
for i=1:N-1
    p=i;
    while Aug(p,p)==0 & p<=N
        p=p+1;
    end
    if p==N
       print('no unique solution exists')
       break
    end
    for j=i+1:N
        m=Aug(j,i)/Aug(i,i);
        Aug(j,:)=Aug(j,:)-m*Aug(i,:);
    end
end
Aug
X(N)=Aug(N,N+1)/Aug(N,N);
for k=N-1:1
    X(k)=(Aug(k,N+1)-Aug(k,k+1:N)*X(k+1:N))/Aug(k,k);
end
X
    
    
