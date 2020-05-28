
K=load('avsLyapunov-less-combination.dat');

amplitude=K(:,2);

L1=K(:,3);

L2=K(:,4);

L1array=zeros(30,1);

L2array=zeros(30,1); 

for block=0:29

	for elements=1:999
	
	amparray(block+1)=mean(amplitude(elements+block*999));
	
	L1array(block+1)=mean(L1(elements+block*999));

	L2array(block+1)=mean(L2(elements+block*999));

	total(block+1)=mean(L1(elements+block*999)+L2(elements+block*999));

	end 
end

f=figure(1);

subplot(3,1,1),plot(amparray,L1array,'.-','color', 'blue'); 

ylabel('$\lambda_1$','interpreter','latex','Fontsize',14),grid on;

% title(['$ latex $'],'interpreter','latex')

xlabel('$a$','interpreter','latex','Fontsize',14);

title('(a)','Fontsize',14);


subplot(3,1,2),plot(amparray,L2array,'.-','color', 'red');

ylabel('$\lambda_2$','interpreter','latex','Fontsize',14),grid on;

xlabel('$a$','interpreter','latex','Fontsize',14);

title('(b)','Fontsize',14);

%subplot(3,1,3),plot(amparray,total,'linewidth',2,'color','green','marker','none');

subplot(3,1,3),plot(amparray,total,'.-','color','m');


ylabel('$\lambda_1$+$\lambda_2$','interpreter','latex','Fontsize',14);grid on;

xlabel('$a$','interpreter','latex','Fontsize',14);

title('(c)','Fontsize',14 );

print -depsc a-vs-lambda-less

%saveas(f,'a-vs-lambda-less.png');grid on;clear all;
