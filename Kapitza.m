clear all;close all;clc;

LL=load('Lyapunov-threshold.dat');

PP=load('Poincare-threshold.dat');

PS=load('Phase-space.dat');


xx=LL(:,2);

%Theta=abs(mod(xx,2*pi));

Theta=mod(xx,2*pi)

Thetadot=LL(:,3);

tt=LL(:,1);

%($2+pi - (floor(($2+pi)/(2*pi))*(2*pi) )):3 w d

theta1=(PS(:,2)+pi-(floor(PS(:,2)/(2*pi))*(2*pi)));

%theta1=mod(PS(:,2),pi);
thetadot1=PS(:,3);

ptheta=PP(:,3);

pthetadot=PP(:,4);

f=figure(1);

T=tt(2)-tt(1);

%Computing Fourier Spectrum

Fs=1/T;       % sampling rate

L=length(tt);

%t=(0:10^5-1)*T;% sampling instants

NFFT = 2^nextpow2(L);

f1 = 2*pi*Fs/2*linspace(0,1,NFFT/2+1);

f11=f1(find(f1<100));

fX=fft(Thetadot-mean(Thetadot),NFFT)/L;        %Complex Amplitudes

fX=abs(fX);

f=figure(1);


%Time series for theta

subplot(3,2,1),plot(tt,Theta,'linewidth', 2, 'color', 'blue', 'marker', 'none'); 

ylabel('$\theta$','interpreter','latex','Fontsize',14);

xlabel('Time','interpreter','latex','Fontsize',10);

%xlim([0 10]);

title('(a)','interpreter','latex','Fontsize',10);

%Timeseries for thetadot

subplot(3,2,2),plot(tt,Thetadot,'linewidth', 1, 'color', 'blue', 'marker','none');

ylabel('$\dot{\theta}$','interpreter','latex','Fontsize',14);

xlabel('Time','interpreter','latex','Fontsize',10);

%xlim([0 10]);

title('(b)','interpreter','latex','Fontsize',10);


%Autocorrelation

auto=autocorr(Theta,5000);

subplot(3,2,3),plot(auto,'linewidth', 1, 'color', 'blue', 'marker','none');

ylabel('Autocorrelation','interpreter','latex','Fontsize',10);

xlabel('Time','interpreter','latex','Fontsize',10);

%set(gca, 'XTick', [0.0 1000 2000 3000],'XTickLabel',{'0' '2*10^3' '3*10**3' '4*10**3' '5*10^3'},'Fontname','symbol');

title('(c)','interpreter','latex','Fontsize',10);

%Fourier Spectrum

subplot(3,2,4),plot(f11,(fX(1:length(f11))))

ylabel('$F(\omega)$','interpreter','latex','Fontsize',10);

xlabel('$\omega$','interpreter','latex','Fontsize',10);

%set(gca, 'XTick', [0.0 1000 2000 3000],'XTickLabel',{'0' '2*10^3' '3*10^3' '4*10^3' '5*10^3'},'Fontname','symbol');

title('(d)','interpreter','latex','Fontsize',10);

%Poincare-Section

subplot(3,2,5),plot(theta1,thetadot1,'.');

ylabel('$\dot{\theta}$','interpreter','latex','Fontsize',14);

xlabel('$\theta$','interpreter','latex','Fontsize',14);

axis tight;



set(gca, 'XTick', [  pi/2 pi  3*pi/2  2*pi    5*pi/2],'XTickLabel',{  'p/2' 'p ' '3p/2' '2p' '5p/2'},'Fontname','symbol');

set(gca, 'YTick', [-3*pi -2*pi -pi 0.0 pi 2*pi 3*pi],'YTickLabel',{'-3p' '-2p' 'p' '0' 'p' '2p' '3p'},'Fontname','symbol');


%xlim([0 4*pi]);

title('(e)','interpreter','latex','Fontsize',10);

%Stroboscopic-plot

subplot(3,2,6),plot(ptheta,pthetadot,'.');

ylabel('$\dot{\theta}$','interpreter','latex','Fontsize',14);

xlabel('$\theta$','interpreter','latex','Fontsize',14);

axis tight;

%xlim([0 2*pi]);

set(gca, 'XTick', [ 0 pi/2  pi  3*pi/2 2*pi],'XTickLabel',{ '0' 'p/2' 'p' '3p/2' '2p' },'Fontname','symbol');

%set(gca, 'XTick', [0.0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0' 'p/2' 'p' '3p/2' '2p'},'Fontname','symbol');

set(gca, 'YTick', [-3*pi -2*pi -pi 0.0 pi 2*pi 3*pi],'YTickLabel',{'-3p' '-2p' 'p' '0' 'p' '2p' '3p'},'Fontname','symbol');

title('(f)','interpreter','latex','Fontsize',10);

%print('image', '-dpng', '-r300');

print -depsc Kapitza-threshold.eps

%saveas(f,'Kapitza-below-threshold.eps');grid on;clear all;
