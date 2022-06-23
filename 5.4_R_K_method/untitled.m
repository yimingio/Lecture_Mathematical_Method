a=0;
b=2;
N=10;
ya=0.5;
f=@(t,y) y-t^2+1;
[t,w]=func_rk4(a,b,ya,N,f);