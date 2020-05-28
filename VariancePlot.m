%variance plot for Kapitza Pendulum (LongTimeSeries)


clear all;close all;clc;


K=load('DifferentCondition.dat');

Theta= K(:,2);

ThetaDot= K(:,4);

T=zeros(100,1);%matrix for number of blocks

TD=zeros(100,1);%matrix for number of blocks

 ThetaVariance=zeros(99999,1);%matrix for elements in each block
 
ThetaDotVariance=zeros(99999,1);%matrix for elements in each block

count=0;

for element=1:99999%number of elements
    
    count=count+1
    
for block=0:99  %number of blocks
    
    T(block+1)=Theta(99999*block+element);
    
    TD(block+1)=ThetaDot(99999*block+element);
    
end
 
 ThetaVariance(element)=var(T);%variance
 
 ThetaDotVariance(element)=var(TD);
end




f1=figure(1);

subplot(2,2,1);

plot(ThetaVariance);

title('Thetaplot');

subplot(2,2,2);

loglog(ThetaVariance);

title('loglogplotTheta');

subplot(2,2,3);

plot(ThetaDotVariance);

title('Thetadotplot');

subplot(2,2,4);

loglog(ThetaDotVariance);

title('loglogplotThetadot');

saveas(f1,'variancePlot.jpg')


