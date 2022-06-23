function [t,w] = func_Euler(a,b,ya,N,f)
h=(b-a)/N;
t=[];
w=[];

t(0 +1)=a;
w(0 +1)=ya;
for i=1:1:N
    w(i +1)=w(i-1 +1)+h*f(t(i-1 +1),w(i-1 +1));
    t(i +1)=a+i*h;
end                         %for example w(1)=corressponding to w0
end

