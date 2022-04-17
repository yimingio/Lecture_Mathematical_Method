%Software:      MATLAB R2021b
%Code:          PHY2009481
%Name:          Guo Yiming
%% Q5a Graph Plotting
x=linspace(-10,10,1000);
y=2-x.^2.*sin(x);
plot(x,y,'LineWidth',2,'Color','k');
line([0 0],[-80 80],'color','k');
line([-10 10],[0 0],'color','k');
grid on;