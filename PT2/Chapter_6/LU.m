A=[1,1,0,3;
    2,1,-1,1;
    3,-1,-1,2;
    -1,2,3,-1];

[L,U]=func_LU(A);




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