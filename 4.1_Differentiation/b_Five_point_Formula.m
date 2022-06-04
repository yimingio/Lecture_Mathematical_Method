%% Five point formula


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
funcdiffmid=1/(12*h)*(func(x0-2*h)-8*func(x0-h)+8*func(x0+h)-func(x0+2*h));                                         %Determine the three point end point function
funcdiffend=1/(12*h)*(-25*func(x0)+48*func(x0+h)-36*func(x0+2*h)+16*func(x0+3*h)-3*func(x0+4*h));                   %Determine the three point mid point function
%% Display the result
fprintf("The 5 point end point is %8.4f\n",funcdiffend);
fprintf("The 5 point mid point is %8.4f\n",funcdiffmid);
%% Appendix
function val=func(x)
    val=x*exp(x);
end