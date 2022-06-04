
h=0.05;
x0=1;
fd = 1/(2*h)*(func(x0+h)-func(x0-h));
disp(fd)


function val = func(x)
    val=exp(x);
end

