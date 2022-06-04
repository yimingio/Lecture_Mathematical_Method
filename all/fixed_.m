%% Fixed Point
% 
%% Initialization
clc;
clear;
close all; %Clear before
%% Add initial condition

tol = 1.e-2;            % Tolorence of result
p = 1.5;                  % Initial left value
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
while (itcount <= nmax)      % When the number of iterations does not overflow, the error is still too large
    itcount =itcount+1;
    p = func(p);
    z =[z,p];
    if itcount==nmax
        x_final=p;
    end
end
%% 
% Plot the Graph 绘制图像
% plot(z(1:itcount),zeros(itcount,1),'rx');
%% Show the result
%
fprintf('Converged solution after %d iterations',itcount);
fprintf(' is %.7e\n',x_final);

%% Appendix
%
function y = func(x)
    y = (10/(4+x))^(1/2);
end
