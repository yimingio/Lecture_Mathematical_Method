a=0;
b=1;
y=@(x) exp(x);
n=4;
val =func_comp_midpoint(a,b,n,y);
disp(val);







function XI = func_comp_midpoint(a,b,n,y)
xsample=linspace(a,b,n+3); %create the sample matrix
h=(b-a)/(n+2);
XI1=0;
for i=0:n/2
    XI1=XI1+y(xsample(2*i +2));
end
XI=2*h*XI1;
end
