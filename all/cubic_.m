


x = [0 1 2 3];

y = [1 1.7 2.8 3.9];

coeff=cubic(x,y);
disp(coeff);


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
    xT = x.';
    xT = xT(1:n-1,1);
    coeff=[xT,coeff];
end