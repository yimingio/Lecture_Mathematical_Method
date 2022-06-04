syms x

f = sqrt(x^2+1);
diff_f = diff(f);

y1 = subs(diff_f,x,1);
y2 = subs(diff_f,x,2);
y3 = subs(diff_f,x,3);


funcdiff=matlabFunction(diff_f);
