

all: rocket.exe

rocket.exe: rocket.f
	gfortran rocket.f -o rocket.exe

plot: altitude.png
	gnuplot -persist rocket.gnuplot


