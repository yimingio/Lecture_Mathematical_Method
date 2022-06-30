a=0;
b=2;
N=10;
ya=0.5;
f=@(t,y) y-t^2+1;



[t,w] = func_mod_euler(a,b,ya,N,f);








function [t,w] = func_mod_euler(a,b,ya,N,f)
h=(b-a)/N;
t=[];
w=[];



t(0 +1)=a;
w(0 +1)=ya;
for i=1:1:N
    t(i +1)=a+i*h;
    fi=f(t(i-1 +1),w(i-1 +1));
    fi2=f(t(i +1),w(i-1 +1) + h*fi);
    T=fi+fi2;
    w(i +1)=w(i-1 +1)+h*T/2;

end
end