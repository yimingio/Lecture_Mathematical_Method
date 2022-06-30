x0=2;               %initial point
h=1;             %step size
nrow=4;             %Input the row need
func=@(x) x*exp(x);
hg=[];
N=[];
hg(1)=h;
for i=1:nrow
    Ntemp=1/(2*h)*(func(x0+h)-func(x0-h));
    N(i,1)=Ntemp;
    if i > 1
        for j=2:i
            Ntemp=N(i,j-1)+(N(i,j-1)-N(i-1,j-1))/(4^(j-1)-1);
            N(i,j)=Ntemp;
        end
    end
end

