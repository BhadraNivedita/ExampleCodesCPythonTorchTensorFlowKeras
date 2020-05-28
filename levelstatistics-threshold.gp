set terminal postscript color eps enhanced 'Helvetica, 24'
set output 'levelstat-threshold.eps'
#unset key 
set xlabel 'S( energy level spacing)'
set ylabel 'P(S)'
#set label ' Poissonian distribution fit'
f(x)=a*exp(-b*x*x)+c*x*exp(-d*x**2)
fit f(x)  'Histo_threshold.dat' via a,b,c,d

pl  'Histo_threshold.dat' w l linewidth 2  ti 'raw data',f(x) w l linewidth 2 lt -1  ti 'f(x)=a exp(-bx^2)+c x exp(-dx^2)'
