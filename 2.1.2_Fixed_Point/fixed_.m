%% Fixed Point
% 
%% Initialization
clc;
clear;
close all; %Clear before
%% Add initial condition

tol = 1.e-3;            % Tolorence of result
p = 2;                  % Initial left value
nmax = 100;             % Number of initialization iterations
itcount = 0;            % Number of initialization iterations
error = 10000;          % Error initialization 
z=[];
%% 
% Plot the graph

% xval = linspace(a,b,100);
% for i=1:100
%     fval(i) = func(xval(i));        % Calculate function image points
% end
% plot(xval,fval);                    % Plot the function
% z
% grid on;
% hold on;
%% Iteration
% 
while (itcount <= nmax && error >=tol)      % When the number of iterations does not overflow, the error is still too large
    itcount =itcount+1;
    p = func(p);
    z =[z,p];
    if (error < tol)
        p_final = p;
    end
end
%% 
% Plot the Graph 绘制图像
plot(z(1:itcount),zeros(itcount,1),'rx');
%% Show the result
%
if (itcount <nmax)
    val = func(x);
    fprintf('Converged solution after %d iterations',itcount);
    fprintf(' is %.7e\n',x_final);
else
    fprinf('Not converged after %5d iterations',nmax);
end

%% Appendix
%
function y=func(x)
    y = exp(x)-9.97418;
end