function [midresult,endresult]=func_3pfdiff(x0,h,func)
    endresult=1/(2*h)*(-3*func(x0)+4*func(x0+h)-func(x0+2*h));    %Determine the three point end point function
    midresult=1/(2*h)*(func(x0+h)-func(x0-h));                    %Determine the three point mid point function
end