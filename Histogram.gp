set terminal postscript  eps enhanced color font 'Helvetica,10'

f(x)=x*exp(-(x+1.36)**2)
set yrange[0.00:0.08]
set xrange[0.0:5]
pl 'Histo_above.dat' w lp linewidth 5
set output 'Levelstatistics.eps'
repl f(x) w lp linewidth 4
