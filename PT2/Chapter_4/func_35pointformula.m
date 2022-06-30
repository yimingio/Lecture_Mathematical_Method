%% Initialization
close all;clear;clc;
%% Input data
x0=2;         %initial point
h= 0.1;         %step size
func=@(x) x*exp(x);
num_points=3;
type="mid";
result=f35pf(x0,h,func,num_points,type);
disp(result)

%% Function
function result=f35pf(x0,h,func,num_points,type)
if num_points == 3
    if type == "mid"
        funcdiffmid=1/(2*h)*(func(x0+h)-func(x0-h));                    %Determine the three point mid point function
        result=funcdiffmid;
    elseif type == "end"
        funcdiffend=1/(2*h)*(-3*func(x0)+4*func(x0+h)-func(x0+2*h));    %Determine the three point end point function
        result=funcdiffend;
    else
        fprintf("Error input the type name please input 'mid' or 'end'\n");
    end
elseif num_points ==5
    if type == "mid"
        funcdiffmid=1/(12*h)*(func(x0-2*h)-8*func(x0-h)+ ...
            8*func(x0+h)-func(x0+2*h));                                         %Determine the three point end point function                   %Determine the three point mid point function
        result=funcdiffmid;
    elseif type == "end"
        funcdiffend=1/(12*h)*(-25*func(x0)+48*func(x0+h)- ...
            36*func(x0+2*h)+16*func(x0+3*h)-3*func(x0+4*h));                   %Determine the three point mid point function
        result=funcdiffend;
    else
        fprintf("Error input the type name please input 'mid' or 'end'\n");
    end
else
    fprintf("Error please input 3 or 5\n");
end
end