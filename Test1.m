N =50;

w=40;

P=0.001;

g=9.81;

a=0.11;

phi=pi;

in=1:N;

%for phi=-pi/2.0 :.5:pi/2.0
%fileID = fopen('Histo_above-flux-piby2.dat','w');
%H0 = eye(N);
%H0( 1 : ( N+1 ) : end ) = in*(w)+(1.0/h)*((P*P)/2.0-g*cos(phi+theta)-(1.0/4.0)*a*a*w*w*cos(phi+theta)*cos(phi+theta));
H0=zeros(N);

%H0 = H0 + eye(N);%Diagonal elements

    H0( 2 : ( N+1 ) : end ) = -0.5*a*g*ones(1,N-1);
    H0( N+1 : ( N+1 ) : end ) = -0.5*a*g*ones(1,N-1);
    
    %H0(1,N) =(1.0/h)*(-1i*a*w*P*sin(phi+theta)/2.0-(a*g)/2.0);
    %H0(N,1) = (1.0/h)*(+1i*a*w*P*sin(phi+theta)/2.0-(a*g)/2.0);
    
    
    H0( 3 : ( N+1 ) : end-N ) = 0.125*a*a*w*w*cos(phi)*cos(phi)*ones(1,N-2) ;
    H0( 2*N+1 : ( N+1 ) : end ) = 0.125*a*a*w*w*cos(phi)*cos(phi)*ones(1,N-2);
    
    
    %H0(1,N-1)=(1.0/h)*(1.0/8.0)*a*a*w*w*cos(phi+theta)*cos(phi+theta);
    %H0(N-1,1)=(1.0/h)*(1.0/8.0)*a*a*w*w*cos(phi+theta)*cos(phi+theta);
    %H0(2,N)=(1.0/h)*(1.0/8.0)*a*a*w*w*cos(phi+theta)*cos(phi+theta);
    %H0(N,2)=(1.0/h)*(1.0/8.0)*a*a*w*w*cos(phi+theta)*cos(phi+theta);
    
  %subplot(211);
  
  plot(eig(H0),'.')%eigenvalue vs pi plot
 % eig(H0)
  %min_eigen_value=min(eig(H0));
  %hold on;
%end

%dEv = diff(sort(eig(H0)));
%[x,y]=hist(dEv,50);

%xx=x/sum(x);
%yy=y/sum(y.*xx);
%for ii=1:length(xx)

%save('pqfile.txt','xx','yy','-ascii')
%fprintf(fileID,'%2.6f %2.6f\n',yy(ii),  xx(ii));
%subplot(212);
%end 
%plot(yy,xx,'o-')
%fclose(fileID);
%print -depsc -painters LevelStatistics(Histogram_omega40_above_flux_piby2).eps;


%A  = rand(10,1);
%B = rand(10,1);
%header1 = 'Hello';
%header2 = 'World!';
%fid=fopen('MyFile.txt','w');
%fprintf(fid, [ header1 ' ' header2 '\n']);
%fprintf(fid, '%f %f \n', [A B]');
%fclose(fid);
