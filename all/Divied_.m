%% Newton’s Divided-Difference Formula
% 
%% 
% Clear Previous Data

clc;
clear;
%% 
% Input the initial point:

X = [1.0 1.3 1.6 1.9 2.2]; % Initial Value
n = size(X);
n = n(2);
F=[];
Y=[0.7651977 0.6200860 0.4554022 0.2818186 0.1103623];
% for i = 1:n
%     Y(i)=func(X(i));   %Exact value for Y
% end
%% Iteration Rule
% Start build the interated Interpolation:
% 
% $$F_{i, j}=\frac{F_{i, j-1}-F_{i-1, j-1}}{x_{i}-x_{i-j}}$$

for i = 1:n
    F(i,1)=Y(i); %Input the value in Q(i+1,0+1)
end

for i=2:n
    for j=2:i
        F(i,j)=(F(i,j-1)-F(i-1,j-1))/(X(i)-X(i-j+1));
    end
end
%% Appendix Function

% function y = func(x)
%     y=x^2;
% end