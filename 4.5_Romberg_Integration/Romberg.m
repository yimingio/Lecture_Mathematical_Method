a=0;
b=4;
f=@(x) x/(sqrt(x^2+9));
n=5;

Rresult = func_romberg(a,b,n,f)




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
