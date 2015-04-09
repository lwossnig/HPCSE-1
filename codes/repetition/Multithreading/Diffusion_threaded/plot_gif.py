out = open("plot.gnu", "w\n")
out.write("reset\n")
out.write("set term gif animate\n")
out.write("set output \"animate.gif\"\n")
out.write("set style data dots\n")
out.write("set xrange [-1.0:1.0]\n")
out.write("set yrange [-1.0:1.0]\n")
out.write("set zrange [-1.0:1.0]\n")
out.write("set hidden3d\n")
out.write("set grid\n")
out.write("set size square\n")
for i in range(1000):
    out.write("splot \"density_"+str(i)+".dat\" u 1:2:3 w l\n")