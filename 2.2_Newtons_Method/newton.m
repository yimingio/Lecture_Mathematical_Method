%% Newton Method
% 
%% Initialization

clc;
clear;
close all; % Clear before
%% Add initial condition

tol = 1.e-10;           % Tolorence of result
nmax = 100;             % Maximum number of iterations
p_1 = 3;                % Initial Prediction value
xmin = 0;               % Min x of plot interval
xmax = pi/2;            % Max x of plot interval
itcount = 0;            % Number of initialization iterations
error = 10000;          % Error initialization
%% Plot the graph
% About the drawing related parameters, and start the graphic diagram to draw 
% simple pictures


xval = linspace(xmin,xmax,100); 
for i=1:100
    fval(i) = func(xval(i));        % Calculate function image points
end
plot(xval,fval);                    % Plot the function
%% 
% Start point calculation:

grid on;
hold on;
while (itcount <= nmax && error >=tol)  % When the number of iterations does not overflow, the error is still too large
    itcount =itcount+1;
    p = p_1-func(p_1)/funcdiff(p_1);    % Iteration function of the Newton's Method
    z(itcount) = p;                     % Record the result into the dataset
    error = abs(p-p_1);                 % Calculate the error
    if (error < tol)                    % Campare the result 
        p_final = p;                    % If the result is acceptable then record it as the final result
    else
        p_1=p;                          % Set the result as the initial value for nex iteration
    end
end
%% 
% Plot the Graph 绘制图像

plot(z(1:itcount),zeros(itcount,1),'rx'); % Plot the point with the dataset
if (itcount <nmax)
    val = func(p);
    fprintf('Converged solution after %d iterations',itcount);
    fprintf(' is %.7e\n',p_final);
else
    fprinf('Not converged after %5d iterations',namx);
end
%% Appendix: Function
% Here is an Example of the function

function val = func(x)
    val = cos(x)-x;
end
%% 
% Here is the derivative of the above function

function val = funcdiff(x)
    val = -sin(x)-1;
end