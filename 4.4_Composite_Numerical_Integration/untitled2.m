clc;close all;clear;
a=0;
b=2;
N=10;
h= (b-a)/N;
t(1)=a;
y(1)=0.5;
dy=@(t, y) y-t^2+1;
for i=1:N
    t(i+1)=t(i)+h;
    y(i+1)=y(i)+h*dy(t(i),y(i));
end
ty=[t;y];
fprintf('%f   %f \n ',ty);