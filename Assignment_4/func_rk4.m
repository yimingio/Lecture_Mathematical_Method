function [t,w,k] = func_rk4(a,b,ya,N,f)
h=(b-a)/N;
t=[];
w=[];
k1=[];
k2=[];
k3=[];
k4=[];


t(0 +1)=a;
w(0 +1)=ya;
for i=1:1:N
    t(i +1)=a+i*h;

    k1d=h*f(t(i-1 +1),w(i-1 +1));
    k2d=h*f(t(i-1 +1)+h/2,w(i-1 +1)+1/2*k1d);
    k3d=h*f(t(i-1 +1)+h/2,w(i-1 +1)+1/2*k2d);
    k4d=h*f(t(i-1 +1) +h,w(i-1 +1)+k3d);
    k1(i +1)=k1d;
    k2(i +1)=k2d;
    k3(i +1)=k3d;
    k4(i +1)=k4d;
    w(i +1)=w(i-1 +1)+1/6*(k1d+2*k2d+2*k3d+k4d);
    
end
k=[k1;k2;k3;k4];
k=k';
