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