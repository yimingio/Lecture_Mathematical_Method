
a=0;
b=4;
y=@(x) x/sqrt(x^2+9);
n=16;
val =func_comp_trapezoidal(a,b,n,y);
disp(val);




function val = func_comp_trapezoidal(a,b,n,y)
xsample=linspace(a,b,n+1); %create the sample matrix
h=(b-a)/n;
XI0=y(a)+y(b);
XI1=0;
for i=1:n-1
    XI1=XI1+y(xsample(i+1));
end
val=h*(XI0+2*XI1)/2;
end