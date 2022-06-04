clc;clear;close all;

a=1;
b=4;
y=@(x) sqrt(x);
nromberg=5;
nsimpson=2^nromberg;

simpsonresult=func_comp_simpson(a,b,nsimpson,y);
rombergresult=func_romberg(a,b,nromberg,y);

simpsonerror=abs(simpsonresult-14/3);
rombergerror=abs(rombergresult(nromberg,nromberg)-14/3);
simpsonrerror=simpsonerror/(14/3);
rombergrerror=rombergerror/(14/3);

fprintf("The result obtained by Composite Simpson Formula is %8.8f,\n and the absolute error is %8.8f,\n relative error is %8.8f\n\n\n\n",simpsonresult,simpsonerror,simpsonrerror);
fprintf("The matrix obtained by Romberg integration is\n");
disp(rombergresult);

fprintf("and the absolute error is %8.8f,\n relative error is %8.8f\n\n",rombergerror,rombergrerror);


% Do some comment
simpsonerrorg=[];
rombergerrorg=[];
for i=1:10
    nromberg=i;
    nsimpson=2^nromberg;
    
    simpsonresult=func_comp_simpson(a,b,nsimpson,y);
    rombergresult=func_romberg(a,b,nromberg,y);
    
    simpsonerror=abs(simpsonresult-14/3);
    rombergerror=abs(rombergresult(nromberg,nromberg)-14/3);
    simpsonrerror=simpsonerror/(14/3);
    rombergrerror=rombergerror/(14/3);
    simpsonerrorg=[simpsonerrorg,simpsonrerror];
    rombergerrorg=[rombergerrorg,rombergrerror];
end
x=linspace(1,10,10);
plot(x,simpsonerrorg,x,rombergerrorg,LineWidth=2);
legend("Composite Simpson Formula","Romberg Integration")
xlim([1,10]);
ylim([-0.005,0.04]);
title("The relative error for Composite Simpson Formula and Romberg Formula of $\int_1^4 \sqrt{t} \mathrm d t$ with $2^n$ interval",Interpreter="latex",FontSize=20);

























function val=func_comp_simpson(a,b,n,y)
    xsample=linspace(a,b,n+1); %create the sample matrix
    h=(b-a)/n;
    XI0=y(a)+y(b);
    XI1=0;
    XI2=0;
    for i=1:n-1
        if mod(i,2) == 0
            XI2=XI2+y(xsample(i+1));
        else
            XI1=XI1+y(xsample(i+1));
        end
    end
    val=h*(XI0+2*XI2+4*XI1)/3;
end


function Rresult = func_romberg(a,b,n,f)
    h=b-a;
    R(1,1)=h/2*(f(a)+f(b));
    Rresult(1,1)=R(1,1);
    for i=2:n
        sum=0;
        for k=1:2^(i-2)
            sum=sum+f(a+(k-0.5)*h);
        end
        R(2,1)=1/2*(R(1,1)+h*sum);
        for j=2:i
            R(2,j)=R(2,j-1)+(R(2,j-1)-R(1,j-1))/(4^(j-1)-1);
        end
        
        h=h/2;
        for j=1:i
            R(1,j)=R(2,j);
            Rresult(i,j)=R(2,j);
        end
    end
end


