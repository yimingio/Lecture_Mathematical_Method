%% Matlab
% PHY2009481 Q1


%% Bisection Method
% 
%% Initialization
%
clc;
clear;
close all; %Clear before
%% Add initial condition
%
tol = 5.e-9;            % Tolorence of result
a = 0;                  % Initial left value
b = 1;                  % Initial right value
nmax = 100;             % Number of initialization iterations
itcount = 0;            % Number of initialization iterations
error = b-a;          % Error initialization    
ag=[];
bg=[];
fag=[];
fxg=[];
xg=[];
errorg=[];
%% Plot the Function graph
% 
% xval = linspace(a,b,100);
% for i=1:100
%     fval(i) = func(xval(i));        % Calculate function image points
% end
% plot(xval,fval);                    % Plot the function
%grid on;
% hold on;
%% Iteration(!!!!This part no need to change)
% 
while (itcount <= nmax && error >=tol)      % When the number of iterations does not overflow, the error is still too large
    ag=[ag;a];
    bg=[bg;b];
    itcount =itcount+1;
    x=(a+b)/2;                              % Iteration function of the Bisection Method
    xg=[xg;x];
    z(itcount)=x;                           % Record the result into the dataset
    fa=func(a);
    fx=func(x); 
    fb=func(b);
    fag=[fag;fa];
    fxg=[fxg;fx];
    error=(1-0)/(2^itcount);                %From Sauer Book
    errorg=[errorg;error];
    if (error < tol)
        x_final = x;
    elseif (fa*fx <0)
        b = x;
    else
        a = x;
    end
end
ig=linspace(1,itcount,itcount);
ig=ig.';
%% Show the result
result=[ig,ag,bg,xg,fag,fxg,errorg];

format shortG

fprintf("Iteration times i=%d\n",itcount);
fprintf("Next is the list for each number\n");
fprintf("            i           ai           bi           pi        f(ai)         f(pi)        error\n");
disp(result)
%% Plot the Dot Graph 绘制图像
% 
% plot(z(1:itcount),zeros(itcount,1),'rx');
%% Appendix
%
function y = func(x)
    y = exp(x)-x^2+3*x-2;
end


