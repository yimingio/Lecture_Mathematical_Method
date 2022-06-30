%Composite Midpoint rule
%matlab

a=0;
b=1;
f=@(x) exp(x);
n=4;

h=(b-a)/(n+2);
XI1=0;

for i=0:n
    X=a+(i+1)*h;
    if mod(i,2)==0
        XI1=XI1+f(X);
    end
end

XI=2*h*XI1;
fprintf('the result is %11.9f\n',XI)