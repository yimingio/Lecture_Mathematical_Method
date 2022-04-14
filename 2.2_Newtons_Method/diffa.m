function val = diffa(x,tol)
    range = tol/100;
    fx = func(x);
    fx1 = func(x+range);
    val = (fx1-fx)/range;
end

 