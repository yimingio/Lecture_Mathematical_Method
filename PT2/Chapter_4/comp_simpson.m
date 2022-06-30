a=0;
b=1;
y=@(x) x^2;
n=4;
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