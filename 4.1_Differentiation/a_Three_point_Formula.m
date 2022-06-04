%% Three point formula


%% Initialization
close all;clear;clc;
%% Input data
x0=2;         %initial point
h= 0.1;         %step size
% Function is written in the Appendix
%% Exact function Generator
% syms x
% f=x*exp(x);
% f_diff=diff(f);
% funcdiff=matlabFunction(f_diff);
% f_2diff=diff(f_diff);
% func2diff=matlabFunction(f_2diff);
%% Processing the data
funcdiffend=1/(2*h)*(-3*func(x0)+4*func(x0+h)-func(x0+2*h));    %Determine the three point end point function
funcdiffmid=1/(2*h)*(func(x0+h)-func(x0-h));                    %Determine the three point mid point function
%% Display the result
fprintf("The 3 point end point is %8.4f\n",funcdiffend);
fprintf("The 3 point mid point is %8.4f\n",funcdiffmid);
%% Appendix
function val=func(x)
    val=x*exp(x);
end