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
        fprintf("Here find tha maximum on $E_{%d}$",p);
        fprintf(". Switch the row  $(E_{%d}) \\leftrightarrow (E_{%d})$ \n",p,i);
    end

    %Generate the matrix we need
    Aaucopy=[];

    for k=1:n
        Aaucopy=[Aaucopy;Aau(NROW(k),:)];
    end
    display_aug_mat(Aaucopy);
    fprintf("Perform the operation \n")
    for j=i+1:n
        m(NROW(j),i)=Aau(NROW(j),i)/Aau(NROW(i),i);
        m(NROW(j),i)=round(m(NROW(j),i),2);
        Aau(NROW(j),:)=Aau(NROW(j),:)-m(NROW(j),i).*Aau(NROW(i),:);
        Aau=round(Aau,2);
        fprintf(" $ ( E_{%d} - (%4.2f) E_{%d}) \\rightarrow (E_{%d}) $, ",NROW(j),m(NROW(j),i),NROW(i),NROW(j));
    end
    %Generate the matrix we need
    Aaucopy=[];

    for k=1:n
        Aaucopy=[Aaucopy;Aau(NROW(k),:)];
    end
    Aaucopy=round(Aaucopy,2);
    fprintf("\n \n");
    display_aug_mat(Aaucopy);
    Ag{i +1}=Aau;
    
 
end


if Aau(NROW(n),n)== 0
    disp("No unique Solution Exist");
else

    fprintf("Finally, the matrix is converted back into a linear system " + ...
        "that has a solution equivalent to the solution of the original system " + ...
        "and the backward substitution is applied: \n \n")
    % Start the backward substituition
    x=[];
    x(n)=Aau(NROW(n),n+1)/Aau(NROW(n),n);
        xn=round(x(n),2);
    fprintf("$$ \n x_{%d}=%4.2f \n $$ \n",n,xn)
    for i=n-1:-1:1
        sum=0;
        for j=i+1:n
            sum=sum+Aau(NROW(i),j)*x(j);
        end
        x(i)=(Aau(NROW(i),n+1)-sum)/Aau(NROW(i),i);
        x(i)=round(x(i),2);
        fprintf("$$ \n x_{%d}=\\left[a_{%d, %d}-\\sum_{j=%d}^{%d} a_{%d j} x_{j}\\right] / a_{%d %d}=%4.2f \n $$ \n",i,i,n+1,i+1,n,i,i,i,x(i));
    end
end




