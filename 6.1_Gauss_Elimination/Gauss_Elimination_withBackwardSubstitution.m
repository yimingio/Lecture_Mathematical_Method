%Gauss Elimination with backward substitution
%matlab

clc;clear;close all;
A=[1,-1,2,-1;
    2,-2,3,-3;
    1,1,1,0;
    1,-1,4,3];
B=[-8;-20;-2;4];
Aug=[A B];
N=length(B);
X=zeros(N,1);
for i=1:N-1
    p=i;
    while Aug(p,i)==0 & p<=N+1
        p=p+1;
    end
    if p==N+1
       print('no unique solution exists')
       break
    end
    Atemp=Aug(i,:);
    Aug(i,:)=Aug(p,:);
    Aug(p,:)=Atemp;
    for j=i+1:N
        m(j,i)=Aug(j,i)/Aug(i,i);
        Aug(j,:)=Aug(j,:)-m(j,i).*Aug(i,:);
    end
    Ag{i}=Aug;
end

if Aug(N,N)== 0
    disp("No unique Solution Exist");
else
end    
X(N)=Aug(N,N+1)/Aug(N,N);
for k=N-1:1
    X(k)=(Aug(k,N+1)-Aug(k,k+1:N)*X(k+1:N))/Aug(k,k);
end
X
    
    
