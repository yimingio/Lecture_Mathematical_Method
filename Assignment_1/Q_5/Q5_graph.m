x=linspace(-10,10,1000);
y=2-x.^2.*sin(x);
plot(x,y,'LineWidth',1);
grid on;
title("$f(x)=2-x^2 \sin(x)$","Interpreter","latex","FontSize",24)