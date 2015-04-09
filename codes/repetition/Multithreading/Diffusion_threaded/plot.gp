set term png
set hidden3d
set grid
set output "before.png"
splot "density_0.dat" u 1:2:3 w l
set output "after.png"
splot "density_serial.dat" u 1:2:3 w l

