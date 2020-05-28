  % Find the eigenvalue and sort them and plot it

a=0.2; g=9.81; w=40.0;

t1=0.5*a*g;t11=1.0-0.4;phi =pi;

t2=0.125*a^2*w^2*cos(phi)*cos(phi);t22=0.4;

k=[-1:0.001:1]

Energy=-2.0*(t1*cos(pi*k)+t2*cos(2*pi*k));

Energy1=-2*(t11*cos(pi*k)+t22*cos(2*pi*k));

%plot(diff(sort(Energy)),'k');

%hold on 

plot(diff(sort(Energy1)),'cyan');

print -depsc -painters energyeigenvaluet1-t2.eps;