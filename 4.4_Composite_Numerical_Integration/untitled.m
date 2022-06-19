%Composite Midpoint rule
%matlab

clc;clear;close all;
a=0;         %Left side point
b= pi;         %Right side point
n= 10;        %subinterval, even
f=@(x) sin(x);

h=(b-a)/(n+2);
XI1=0;

for i=0:n
    X=a+(i+1)*h;
    if mod(i,2)==0
        XI1=XI1+f(X);
    end
end

XI=2*h*XI1;
fprintf('the result is %11.9f\n',XI)
