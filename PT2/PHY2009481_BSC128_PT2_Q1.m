%MatLab R2022a
%PHY2009481 Q1
clc;clear;close all;
a=0;
b=2;
ya=0.5;
n1=2/0.2;
n2=2/0.1;

f=@(t,y) y-t^2+1;
y=@(t) (t+1)^2-0.5*exp(t);


[t1,w1]=func_rk4(a,b,ya,n1,f);
[t2,w2]=func_rk4(a,b,ya,n2,f);

%print the (a)
y1=[];
abserror1=[];
fprintf("(a)---------------------------------------\n")
for i=1:n1+1
    y1(i)=y(t1(i));
    abserror1(i)=abs(y1(i)-w1(i));
end
fprintf("i\tt_i\ty_i\t\tw_i\t\t|y_i-w_i|\n");
for i=1:n1+1
    fprintf("%d     %6.2f\t%6.10f\t%6.10f\t%e\n",i-1,t1(i),y1(i),w1(i),abserror1(i));
end

fprintf("(b)---------------------------------------\n")
y2=[];
abserror2=[];
for i=1:n2+1
    y2(i)=y(t2(i));
    abserror2(i)=abs(y2(i)-w2(i));
end
fprintf("i\tt_i\ty_i\t\tw_i\t\t|y_i-w_i|\n");
for i=1:n2+1
    fprintf("%d     %6.2f\t%6.10f\t%6.10f\t%e\n",i-1,t2(i),y2(i),w2(i),abserror2(i));
end



function [t,w] = func_rk4(a,b,ya,N,f)
h=(b-a)/N;
t=[];
w=[];



t(0 +1)=a;
w(0 +1)=ya;
for i=1:1:N
    t(i +1)=a+i*h;

    k1=h*f(t(i-1 +1),w(i-1 +1));
    k2=h*f(t(i-1 +1)+h/2,w(i-1 +1)+1/2*k1);
    k3=h*f(t(i-1 +1)+h/2,w(i-1 +1)+1/2*k2);
    k4=h*f(t(i-1 +1) +h,w(i-1 +1)+k3);

    w(i +1)=w(i-1 +1)+1/6*(k1+2*k2+2*k3+k4);
end
end











