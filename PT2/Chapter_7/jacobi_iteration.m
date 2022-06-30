A=[10,-1,2,0;
    -1,11,-1,3;
    2,-1,10,-1;
    0,3,-1,8];
B=[6;25;-11;15];
error=10^(-3); %The error
N=1000;        %The maximum iteration times
x=[0;0;0;0];    %Initial value
[x,xg,errorg] = fun_GI(A,B,x,N,error)





function [x,xg,errorg] = fun_GI(A,B,x,N,error)
errorvalue=1000;
xg={};
errorg={};
xg{0 +1}=x;     % The initial value input to cell
%% Generate the D matrix
D=[];
[n,m]=size(A);
for i=1:n
    D(i,i)=A(i,i);
end
%% Generate the L matrix
L=zeros(n);
for i=2:n
    for j=1:i-1
        L(i,j)=-A(i,j);
    end
end

%% Generate the U matrix
U=zeros(n);
for i=1:n-1
    for j=i+1:n
        U(i,j)=-A(i,j);
    end
end
%% Construct the T matrix
T=inv(D)*(L+U);
c=inv(D)*B;
%% Start the iteration
i=0; % Iteration time initizalization
while i<N & errorvalue>error
    xk=T*x+c;
    i=i+1;
    %errorvalue=norm(xk-x,"inf");    %absolute value
    errorvalue=norm(xk-x,"inf")/norm(xk,"inf"); %relative error
    errorg{i+1}=errorvalue;
    xg{i +1}=xk;
    x=xk;

end

%% Display the result
fprintf("k\tx_1(k)\t\tx_2(k)\t\tx_3(k)\t\tx_4(k)\n");
[a,b]=size(xg);
for i=1:b
    xtemp=xg{i};
    fprintf("%d\t",i-1);
    for j=1:n
    fprintf("%8.4f\t",xtemp(j));
    end
    fprintf("\n")
end
end

