set terminal postscript color eps enhanced
set output 'Poincare-threshold.png'
unset key
#set xlabel 'Theta'
#set ylabel 'Thetadot'

pl 'Poincare.dat' u 3:4 w d
