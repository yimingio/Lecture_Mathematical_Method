A=[3,-1;-1,2];
B=[5;4];
error=10^(-3); %The error
N=1000;        %The maximum iteration times
x=[0;0];    %Initial value

[x,xg,errorg] = fun_GI(A,B,x,N,error)