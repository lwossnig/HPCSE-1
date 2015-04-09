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
splot "density_1.dat" u 1:2:3 w l
splot "density_2.dat" u 1:2:3 w l
splot "density_3.dat" u 1:2:3 w l
splot "density_4.dat" u 1:2:3 w l
splot "density_5.dat" u 1:2:3 w l
splot "density_6.dat" u 1:2:3 w l
splot "density_7.dat" u 1:2:3 w l
splot "density_8.dat" u 1:2:3 w l
splot "density_9.dat" u 1:2:3 w l
splot "density_10.dat" u 1:2:3 w l
splot "density_11.dat" u 1:2:3 w l
splot "density_12.dat" u 1:2:3 w l
splot "density_13.dat" u 1:2:3 w l
splot "density_14.dat" u 1:2:3 w l
splot "density_15.dat" u 1:2:3 w l
splot "density_16.dat" u 1:2:3 w l
splot "density_17.dat" u 1:2:3 w l
splot "density_18.dat" u 1:2:3 w l
splot "density_19.dat" u 1:2:3 w l
