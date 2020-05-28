set terminal png
pl 'Lyapunov.dat' u 4 w l title "Liyapunov exponent1"
set output 'Lyapunovexponent.png'
repl 'Lyapunov.dat' u 5 w l   title "Liyapunov exponent2"

