set terminal postscript color eps enhanced 'Helvetica,24'
set output 'levelstat-above.eps'
#unset key 
set xlabel 'S( energy level spacing)'
set ylabel 'P(S)'
#set label ' Poissonian distribution fit'
f(x)=a*exp(-b*x*x)+c*x*exp(-d*x**2)
fit f(x)  'Histo_above.dat' via a,b,c,d
#set yrange[0:1]
pl  'Histo_above.dat' w l linewidth 2  ti 'raw data',f(x) w l linewidth 2 lt -1 ti 'f(x)=a x exp(-b x^2)+c x exp(-d x^2)'
