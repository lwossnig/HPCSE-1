reset
set term gif animate
set output "animate.gif"
set style data dots
set xrange [-1.0:1.0]
set yrange [-1.0:1.0]
set zrange [-1.0:1.0]
set hidden3d
set grid
set size square
splot "density_0.dat" u 1:2:3 w l
splot "density_500.dat" u 1:2:3 w l
splot "density_1000.dat" u 1:2:3 w l
splot "density_1500.dat" u 1:2:3 w l
splot "density_2000.dat" u 1:2:3 w l
splot "density_2500.dat" u 1:2:3 w l
splot "density_3000.dat" u 1:2:3 w l
splot "density_3500.dat" u 1:2:3 w l
splot "density_4000.dat" u 1:2:3 w l
splot "density_4500.dat" u 1:2:3 w l
splot "density_5000.dat" u 1:2:3 w l
splot "density_5500.dat" u 1:2:3 w l
splot "density_6000.dat" u 1:2:3 w l
splot "density_6500.dat" u 1:2:3 w l
splot "density_7000.dat" u 1:2:3 w l
splot "density_7500.dat" u 1:2:3 w l
splot "density_8000.dat" u 1:2:3 w l
splot "density_8500.dat" u 1:2:3 w l
splot "density_9000.dat" u 1:2:3 w l
splot "density_9500.dat" u 1:2:3 w l
