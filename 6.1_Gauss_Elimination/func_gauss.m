function [x,Ag] = func_gauss(A,B)
Ag={};                                                  %Create a cell array to store each matrix
Aau=[A,B];                                              %Generate the augmented matrix
Ag{1}=Aau;                                              %Store the first initial matrix
[n,b]=size(A);
m=[];
for i=1:n-1
    %Check and switch to check pivot to be not 0
    p=i;
    while Aau(p,i)==0 & p <= n-1    
        p=p+1;
    end
    if p==n & Aau(p,i)==0
        disp("No unique solution found, please check");
        break
    end
    Atemp=Aau(i,:);
    Aau(i,:)=Aau(p,:);
    Aau(p,:)=Atemp;
    for j=i+1:n
        m(j,i)=Aau(j,i)/Aau(i,i);
        Aau(j,:)=Aau(j,:)-m(j,i).*Aau(i,:);
        Aau=round(Aau,2);
    end
    Ag{i +1}=Aau;
end

if Aau(n,n)== 0
    disp("No unique Solution Exist");
else



    % Start the backward substituition
    x=[];
    x(n)=Aau(n,n+1)/Aau(n,n);
    x(n)=round(x(n),2);
    for i=n-1:-1:1
        sum=0;
        for j=i+1:n
            sum=sum+Aau(i,j)*x(j);
        end
        x(i)=(Aau(i,n+1)-sum)/Aau(i,i);
        x(i)=round(x(i),2);
    end
end

end

