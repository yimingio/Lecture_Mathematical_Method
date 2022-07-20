%Runge Kutta 4 order method
%matlab

clc;clear;close all;

f = @(t,y) y-t^2+1;    
y=@(t) (t+1)^2-0.5*exp(t);

a = 0;  %initial point
b = 2;  %last point
N1=10; 
N2=20;
y0=0.5;  %initial value


[t1,w1]=RK4(a,b,f,N1,y0);
[t2,w2]=RK4(a,b,f,N2,y0);

%(a)
fprintf("i\tt_i\ty_i\t\tw_i\t\t|y_i-w_i|\n");
for i=1:N1+1
    y1(i)=y(t1(i));
    abs1(i)=abs(y1(i)-w1(i));
    fprintf("%d   %8.2f\t  %6.4f\t  %8.7f\t  %8.7f\n",i-1,t1(i),y1(i),w1(i),abs1(i));
end

%(b)
fprintf("i\tt_i\ty_i\t\tw_i\t\t|y_i-w_i|\n");
for i=1:N2+1
    y2(i)=y(t2(i));
    abs2(i)=abs(y2(i)-w2(i));
    fprintf("%d   %4.2f\t  %8.7f\t  %8.7f\t  %8.7f\n",i-1,t2(i),y2(i),w2(i),abs2(i));
end

function [t,w]=RK4(a,b,f,N,y0)
t=[];
w=[];
t(1)=a;
w(1)=y0;
h=(b-a)/N;
ttemp = a;
wtemp = y0;
for i = 1:1:N 
    k1 = h*f(ttemp,wtemp);
    k2 = h*f(ttemp+h/2, wtemp+k1/2);
    k3 = h*f(ttemp+h/2, wtemp+k2/2);
    k4 = h*f(ttemp+h,wtemp+k3);
    wtemp = wtemp+(k1+2*(k2+k3)+k4)/6;
    ttemp = a+i*h;    
    t(i +1)=ttemp;
    w(i +1)=wtemp;
end
end