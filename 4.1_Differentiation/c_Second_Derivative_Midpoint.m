%% Second Derivative Midpoint Formula

%% Initialization
close all;clear;clc;
%% Input Data
x0=2;         %initial point
h= 0.1;         %step size
% Function is written in the Appendix
syms x
f=x*exp(x);
f_diff=diff(f);
funcdiff=matlabFunction(f_diff);
f_2diff=diff(f_diff);
func2diff=matlabFunction(f_2diff);
%% Processing data
func2diff=1/(h^2)*(func(x0-h)-2*func(x0)+func(x0+h)); %Second Derivative Midpoint Formula
%% Display the data
fprintf("The 2nd dericvative midpoint formula is %8.4f\n",func2diff);
%% Appendix
function val=func(x)
    val=x*exp(x);
end