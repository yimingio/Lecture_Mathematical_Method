a=0;
b=1;
f=@(x) exp(x);
n=4;

F3=compmid(f,a,b,n)




function F3=compmid(f,a,b,n)
I=0;
h=(b-a)/(n+2);
for k=1:n/2+1
    xi=a+(2*k-1)*h;
    I=I+2*h*f(xi);
end
F3=I;
end

