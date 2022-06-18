function [t,w] = func_midpoint(a,b,ya,N,f)
h=(b-a)/N;
t=[];
w=[];
t(0 +1)=a;
w(0 +1)=ya;
for i=1:1:N
    fi=f(t(i-1 +1),w(i-1 +1));
    T=f(t(i-1 +1)+h/2,w(i-1 +1)+h/2*fi);
    w(i +1)=w(i-1 +1)+h*T;
    t(i +1)=a+i*h;
end
plot(t,w);
disp(w(end));
end

