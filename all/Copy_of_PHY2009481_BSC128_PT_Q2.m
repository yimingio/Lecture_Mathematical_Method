%% Matlab
% PHY2009481 Q2

%%
clc;clear;close all;
y1=[];
y2=[];
x1 = linspace(-1,1,5);
x2 = linspace(-1,1,9);
for i=1:5
    y1(i)=func(x1(i));
end


for i=1:9
    y2(i)=func(x2(i));
end
%%
coeff1=cubic(x1,y1);
coeff2=cubic(x2,y2);

%%

fprintf("(b)The coefficient is \n")
format shortG
%fprintf("           -1     -0.48314            0       2.7283")
fprintf("           ai           bi           ci           di\n");
disp(coeff2);

xs=[];
ys=[];

for i=1:8
        xs=linspace(x2(i),x2(i+1),50);
        for j=1:50
            ys(j)=coeff2(i,1)+coeff2(i,2)*(xs(j)-x2(i))+coeff2(i,3)*(xs(j)-x2(i))^2+coeff2(i,4)*(xs(j)-x2(i))^3;
        end
        plot(xs,ys,'r--');
        hold on
end










%%
function coeff=cubic(x,y)
    n=length(x);
    v1=0;
    vn=0;
    A=zeros(n,n); % matrix A is nxn
    r=zeros(n,1);
    for i=1:n-1
        dx(i) = x(i+1)-x(i); 
        dy(i)=y(i+1)-y(i); % define the deltas
    end
    for i=2:n-1
        A(i,i-1:i+1)=[dx(i-1) 2*(dx(i-1)+dx(i)) dx(i)]; % load the A matrix
        r(i)=3*(dy(i)/dx(i) - dy(i-1)/dx(i-1)); % right-hand side
    end
    A(1,1) =1;

    A(n,n) = 1;
    coeff=zeros(n,3);
    coeff(:,2)=A\r;
    for i=1:n-1
        coeff(i,3)=(coeff(i+1,2)-coeff(i,2))/(3*dx(i));
        coeff(i,1)=dy(i)/dx(i)-dx(i)*(2*coeff(i,2)+coeff(i+1,2))/3;
    end
    coeff=coeff(1:n-1,1:3);
    yT = y.';
    yT = yT(1:n-1,1);
    coeff=[yT,coeff];
end
%%

function val=func(x)
    val=1/(1+25*x^2);
end