A=[4.3,6.6,-5.3,6.8;
    2.5,-1.2,+6.6,-2.0;
    5.4,2.2,-2.6,3.5;
    -7.2,5.3,-1.3,4.9];
B=[48.81;
    -30.50;
    45.69;
    -18.15];
Ag={};                                                  %Create a cell array to store each matrix
Aau=[A,B];                                              %Generate the augmented matrix
Ag{1}=Aau;                                              %Store the first initial matrix
[n,b]=size(A);
m=[];
for i=1:n-1
    %Check and switch to check pivot to be not 0
    p=i;
    while Aau(p,p)==0 & p <= n    
        p=p+1;
    end
    if p==n
        disp("No unique solution found, please check");
        break
    end
    Atemp=Aau(i,:);
    Aau(i,:)=Aau(p,:);
    Aau(p,:)=Atemp;
    if i ~= p
        frpintf("switch the row %d with %d \n",i,p);
    end

    for j=i+1:n
        m(j,i)=Aau(j,i)/Aau(i,i);
        Aau(j,:)=Aau(j,:)-m(j,i).*Aau(i,:);
        mr(j,i)=round(m(j,i),2);
        fprintf("Perform the operation E(%d) - (%8.5f) " + ...
            "E(%d) -> E(%d) \n",j,mr(j,i),i,j);
        Aau=round(Aau,2);
    end
    Ag{i +1}=Aau;
    disp(Aau);
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

disp(x);

