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
NROW=[];
for i= 1:n                                              %Initialize row pointer
    NROW(i)=i;
end



disp(Aau);



for i=1:n-1
    fprintf("-----------------------------------------------------------\n")
    %Check and switch to check pivot to be not 0
    p=i;
    for j=i:n
        if abs(Aau(NROW(p),i))<abs(Aau(NROW(j),i))
            p=j;
        end
    end
    if abs(Aau(NROW(p),i))==0
        disp("No unique Solution exist");
        break
    end

    if NROW(i) ~= NROW(p)
        NCOPY=NROW(i);
        NROW(i)=NROW(p);
        NROW(p)=NCOPY;
        fprintf("Here find tha maximum on row %d",p);
        fprintf(". Switch the row between %d and %d \n",p,i);
    end

    %Generate the matrix we need
    Aaucopy=[];

    for k=1:n
        Aaucopy=[Aaucopy;Aau(NROW(k),:)];
    end
    disp(Aaucopy);

    for j=i+1:n
        m(NROW(j),i)=Aau(NROW(j),i)/Aau(NROW(i),i);
        Aau(NROW(j),:)=Aau(NROW(j),:)-m(NROW(j),i).*Aau(NROW(i),:);
        Aau=round(Aau,2);
        fprintf("Perform the operation E(%d) - (%8.5f) " + ...
            "E(%d) -> E(%d) \n",NROW(j),m(NROW(j),i),NROW(i),NROW(j));
    end
    %Generate the matrix we need
    Aaucopy=[];

    for k=1:n
        Aaucopy=[Aaucopy;Aau(NROW(k),:)];
    end
    Aaucopy=round(Aaucopy,2);
    disp(Aaucopy);
    Ag{i +1}=Aau;
    
 
end


if Aau(NROW(n),n)== 0
    disp("No unique Solution Exist");
else
    % Start the backward substituition
    x=[];
    x(n)=Aau(NROW(n),n+1)/Aau(NROW(n),n);
        
    for i=n-1:-1:1
        sum=0;
        for j=i+1:n
            sum=sum+Aau(NROW(i),j)*x(j);
        end
        x(i)=(Aau(NROW(i),n+1)-sum)/Aau(NROW(i),i);
    end
    disp(x);
end




