function [t,w]=func_hi_taylor(a,b,ya,N,f,fd1)
    h=(b-a)/N;
    t=[];
    w=[];
    

    t(0 +1)=a;
    w(0 +1)=ya;
    for i=1:1:N
        T=f(t(i-1 +1),w(i-1 +1))+h/2*fd1(t(i-1 +1),w(i-1 +1)); %calculate the T
        w(i +1)=w(i-1 +1)+h*T;
        t(i +1)=a+i*h;
    end
end