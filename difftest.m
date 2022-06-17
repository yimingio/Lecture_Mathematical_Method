syms x

f = sqrt(x^2+1);
diff_f1 = diff(f,x,1);
diff_f2 = diff(f,x,2);


func=matlabFunction(f);
funcdiff=matlabFunction(diff_f);
