
K=load('avsLyapunov-30combination.dat');

amplitude=K(:,2);

L1=K(:,3);

L2=K(:,4);

L1array=zeros(300,1);

L2array=zeros(300,1); 

for block=0:299

	for elements=1:999
	
	amparray(block+1)=mean(amplitude(elements+block*999));
	
	L1array(block+1)=mean(L1(elements+block*999));

	L2array(block+1)=mean(L2(elements+block*999));

	total(block+1)=mean(L1(elements+block*999)+L2(elements+block*999));

	end 
end

f=figure(1);

subplot(3,1,1),plot(amparray,L1array,'linewidth', 2, 'color', 'blue', 'marker', 'none'); 

ylabel('L1','Fontsize',14),grid on;

xlabel('amplitude','Fontsize',20);

title('L1','Fontsize',14);


subplot(3,1,2),plot(amparray,L2array,'linewidth', 2, 'color', 'red', 'marker','none');

ylabel('L2','Fontsize',14),grid on;

xlabel('amplitude','Fontsize',20);

title('L2','Fontsize',14);

subplot(3,1,3),plot(amparray,total,'linewidth',2,'color','green','marker','none');

ylabel('TotalLyapunovExpoenet','Fontsize',14);grid on;

xlabel('amplitude','Fontsize',20);

title('total','Fontsize',14 );

saveas(f,'LyapunovExponent2.png');grid on;clear all;
