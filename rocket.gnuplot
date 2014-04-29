
set multiplot

set title "Altitude vs Time"
set xlabel "Time (sec)"
set ylabel "Altitude (meters)"
set data style lines
set grid
set nokey

plot 'altitude.dat'
plot 'velocity.dat'
plot 'acceleration.dat'
