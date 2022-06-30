clc; clear; close all; %initialization
%% Input the data

a=0; %input the left endpoint
b=1; %input the right endpoint
ya=1; %input the initial value on y(a)
N=10; %input the number of space
f=@(t,y) 1/y^2;
[t,w]=euler(a,b,ya,N,f);
%% Computation
function [t,w]=euler(a,b,ya,N,f)
h=(b-a)/N;
t=[];
w=[];

t(0 +1)=a;
w(0 +1)=ya;
for i=1:1:N
    w(i +1)=w(i-1 +1)+h*f(t(i-1 +1),w(i-1 +1));
    t(i +1)=a+i*h;
end
end