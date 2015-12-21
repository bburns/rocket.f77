
# gnuplot script file
# see http://www.gnuplot.info/


# set multiplot

set title "Altitude vs Time"
set xlabel "Time (sec)"
set ylabel "Altitude (meters)"

set style data lines
set grid
set nokey

plot 'data\altitude.dat'


set terminal png
set out 'altitude.png'
replot




# plot 'data\velocity.dat'
# plot 'data\acceleration.dat'
