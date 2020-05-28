set terminal postscript color eps enhanced 'Helvetica,24'
set output 'levelstat-below.eps'
#unset key 
set xlabel 'S( energy level spacing)'
set ylabel 'P(S)'
#set label ' Poissonian distribution fit'
f(x)=exp(-a*x*x)
fit f(x)  'Histo_below.dat' via a
set yrange[0:1]
pl  'Histo_below.dat' w l linewidth 4 ti 'raw data',f(x) w l linewidth 4 lt -1 ti 'f(x)=exp(-a x^2)' 
