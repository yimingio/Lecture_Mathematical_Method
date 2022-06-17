function XI = func_comp_midpoint(a,b,n,y)
xsample=linspace(a,b,n+4); %create the equal space points
h=(b-a)/(n+2); %calculate the distance
XI1=0; %initialize the sum

for i=0:n/2
    XI1=XI1+y(xsample(2*i+1)); %add up the all value we want
end

XI=2*h*XI1; %product the distance
end

